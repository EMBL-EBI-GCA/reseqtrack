
package ReseqTrack::HiveProcess::BranchableProcess;

use strict;
use ReseqTrack::Tools::Exception qw(throw);

use base ('Bio::EnsEMBL::Hive::Process');
use ReseqTrack::Tools::FileSystemUtils qw(delete_file);

=head2 fetch_input

    Description : Implements fetch_input() interface method of Bio::EnsEMBL::Hive::Process that is used to read in parameters and load data.
                  Here we have nothing to fetch.

=cut

sub fetch_input {
  my ($self) = @_;
  $self->_process_branch_ids;
  $self->_process_parent_branch_ids;
  $self->_process_child_branch_ids;

  my $root_output_dir = $self->param('root_output_dir');

  my $branch_params_in = $self->param('branch_parameters_in');
  my $universal_branch_params_in = $self->param('universal_branch_parameters_in');
  while (my ($param_key, $arg) = each %$universal_branch_params_in) {
    next if exists $branch_params_in->{$param_key};
    $branch_params_in->{$param_key} = $arg;
  }

  while (my ($param_key, $arg) = each %$branch_params_in) {
    my @param_vals;
    my $arrayref = ref($arg) eq 'ARRAY' ? $arg : [$arg];
    foreach my $array_arg (@$arrayref) {
      my $hashref = ref($array_arg) eq 'HASH' ? $array_arg : {key => $array_arg};
      my $branch_data = $self->_branch_meta_data(key => $hashref->{'key'}, descend => $hashref->{'descend'}, ascend => $hashref->{'ascend'});
      push(@param_vals, map {$_->{'meta_value'}} @$branch_data);
      if ($hashref->{'becomes_inactive'}) {
        $self->_data_to_make_inactive([map {$_->{'branch_meta_data_id'}} @$branch_data]);
        foreach my $file ( map {$_->{'meta_value'}} @$branch_data) {
          next if $file !~ /^$root_output_dir/;
          $self->_files_to_delete($file);
        }
      }
    }
    if (scalar @param_vals) {
      $self->param($param_key, scalar @param_vals == 1 ? $param_vals[0] : \@param_vals);
    }
  }

  $self->_process_output_dir;
  $self->_process_labels;
}

=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
}

=head2 write_output

    Description : Implements write_output() interface method of Bio::EnsEMBL::Hive::Process that is used to deal with job's output after the execution.
                  Dataflows the intermediate results down branch 1, which will be routed into 'intermediate_result' table.

=cut

sub write_output {
    my ($self) = @_;

    $self->_inactivate_data();

    my $branch_params_out = $self->param('branch_parameters_out');
    foreach my $arg (grep {ref($_) ne 'HASH'} values %$branch_params_out) {
      $arg = {key => $arg};
    }

    my $branch_id_hash = $self->{'_branch_id_hash'};
    my @branch_ids = values %$branch_id_hash;

    my %output_data;
    while (my ($param_key, $value) = each %{$self->output_this_branch}) {
      my $db_hash = $branch_params_out->{$param_key};
      my $db_key = $db_hash ? $db_hash->{'key'} : $param_key;
      $output_data{$db_key} = $value;
    }
    $self->_add_branch_data(\@branch_ids, \%output_data);

    my $output_child_branches = $self->output_child_branches;
    my $num_child_branches = scalar @$output_child_branches;
    my $child_branch_ids = $self->_add_branches($num_child_branches, ??SYSTEM??);
    foreach my $i (0..$num_child_branches-1) {
      my %db_data;
      DATA_TYPE:
      while (my ($param_key, $value) = each %{$self->output_child_branches->[$i]}) {
        my $db_hash = $branch_params_out->{$param_key};
        my $db_key = $db_hash ? $db_hash->{'key'} : $param_key;
        $db_data{$db_key} = $value;
      }
      $self->_add_branch_data(??BRANCH_IDs??, \%db_data);
      $self->dataflow_output_id( {'branch_id' => [$child_branch_ids->[$i], ??OTHERS??]}, 2);
    }

    foreach my $i (@{$self->flows_this_branch}) {
      $self->dataflow_output_id( {'branch_id' => \@branch_ids}, $i);
    }

    foreach my $file (@{$self->_files_to_delete}) {
      delete_file($file);
    }
}

# Just in case the child class sets disconnect_when_inactive(1) and then throws an error before changing it back
sub DESTROY {
  my $self = shift;
  return if !$self->{'_data_dbc'};
  $self->data_dbc->disconnect_when_inactive(0);
}

sub _process_branch_ids {
  my ($self) = @_;
  my $param_branch_ids = $self->get_param_array('branch_id');
  throw("no branch id") if !@param_branch_ids;
  my $sql = "SELECT branch_system_id FROM branch WHERE branch_id=?";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my %branch_ids;
  foreach my $branch_id (@$param_branch_ids) {
    $sth->bind_param(1, $branch_id);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    my $row = $sth->fetchrow_arrayref();
    throw("did not recognise $branch_id") if !defined $row;
    $branch_ids{$row->[0]} = $branch_id;
  }
  $self->{'_branch_id_hash'} = \%branch_ids;
}

sub _process_output_dir {
  my ($self) = @_;
  my $root_output_dir = $self->param('root_output_dir');
  my $sub_dirs = $self->get_param_array('output_subdir');
  my $output_dir = join('/', $root_output_dir, @$sub_dirs);
  $self->output_dir($output_dir);
}
sub _process_label {
  my ($self) = @_;
  my $analysis_label = $self->param('analysis_label') || $self->analysis->logic_name;
  my $branch_labels = $self->get_param_array('branch_label');
  my $job_name = join('.', @$branch_labels, $analysis_label);
  $self->job_name($job_name);
}

sub _process_parent_branch_ids {
  my ($self) = @_;
  my %parent_branch_ids;
  my $sql = "SELECT parent_branch_id FROM branch WHERE branch_id=?";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  while (my ($system_id, $branch_id) = each %{$self->{'_branch_id_hash'}}) {
    my @branch_parents;
    my $query_branch_id = $branch_id;
    BRANCH:
    while (1) {
      $sth->bind_param(1, $query_branch_id);
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      my $row = $sth->fetchrow_arrayref;
      last BRANCH if !defined $row || !defined $row->[0];
      $query_branch_id = $row->[0];
      push(@branch_parents, $query_branch_id)
    }
    $parent_branch_ids{$system_id} = \@branch_parents;
  }
  $self->{'_parent_branch_id_hash'} = \%parent_branch_ids;
}

sub _process_child_branch_ids {
  my ($self) = @_;
  my $branch_id_hash = $self->{'_branch_id_hash'};
  my %child_branch_ids;
  my $sql = "SELECT branch_id FROM branch WHERE parent_branch_id=? ORDER BY sibling_index";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  while (my ($system_id, $branch_id) = each %$branch_id_hash) {
    my @ordered_branch_ids = __recursively_query_child_branch($branch_id, $sth);
    foreach my $i (0..$#ordered_branch_ids) {
      $child_branch_ids{$system_id}->{$ordered_branch_ids[$i]} = $i;
    }
  $self->{'_child_branch_id_hash'} = \%parent_branch_ids;
}

sub __recursively_query_child_branch {
  my ($branch_id, $sth) = @_;
  my @child_branch_ids;
  $sth->bind_param(1, $branch_id);
  $sth->execute() or die "could not execute: ".$sth->errstr;
  while (my $row = $sth->fetchrow_arrayref) {
    my $child_branch_id = $row->[0];
    push(@child_branch_ids, $child_branch_id, @{__recursively_query_child_branch($child_branch_id, $sth)});
  }
  return \@child_branch_ids;
}

sub _process_branch_params_in {
  my ($self, $branch_params_in) = @_;
  my $branch_id_hash = $self->{'_branch_id_hash'};
  my $sql1 = "
    SELECT p.process_data_id, p.data_value, count(b2.branch_id) 
    FROM process_data p, branch_data b1, branch_data b2
    WHERE p.process_data_id=b1.process_data_id
    AND p.process_data_id=b2.process_data_id
    AND b1.branch_id=?
    AND p.data_key = ?
    AND p.is_active=1
    GROUP BY p.process_data_id";
  my $sql2 = "
    SELECT b.branch_id, b.branch_system_id
    FROM branch b, branch_data bd
    WHERE b.branch_id=bd.branch_id
    AND bd.process_data_id=?  ";
  my $hive_dbc = $self->data_dbc();
  my $sth1 = $hive_dbc->prepare($sql1) or die "could not prepare $sql1: ".$hive_dbc->errstr;
  my $sth2 = $hive_dbc->prepare($sql2) or die "could not prepare $sql1: ".$hive_dbc->errstr;
  while (my ($param_key, $arg) = each %$branch_params_in) {
    my ($db_key, $ascend, $descend) = ref($arg) eq 'HASH'
                  ? @{$arg}{qw(key ascend descend)}
                  : ($arg, 0, 0);
    my @param_vals;
    while (my ($system_id, $branch_id) = each %$branch_id_hash) {
      $sth1->bind_param(1, $branch_id);
    }
  }
}

sub get_param_array {
  my ($self, $key) = @_;
  my $param = $self->param($key);
  return $param if ref($param) eq 'ARRAY';
  return [] if !defined $param;
  return [$param];
}

sub output_child_branches {
  my ($self, $data_hash, $branch_system) = @_;
  if (scalar keys %args) {
    $self->param('_output_child_branches', $self->output_child_branches);
    push(@{$self->param('_output_child_branches')}, \%args);
  }
  return $self->param('_output_child_branches') || [];
}

sub output_this_branch {
  my ($self, %args) = @_;
  if (scalar keys %args) {
    $self->param('_output_this_branch', \%args);
  }
  return $self->param('_output_this_branch') || {};
}

sub flows_this_branch {
  my ($self, $arg) = @_;
  if (scalar @_ >= 2) {
    my $flows = ref($arg) eq 'ARRAY' ? $arg
              : defined $arg ? [$arg]
              : [];
    $self->{'_flows_this_branch'};
  }
  return $self->{'_flows_this_branch'} || [1];
}

sub _add_process_data {
  my ($self, $key_val_pairs) = @_;
  my $sql = "INSERT INTO process_data (data_key, data_value, is_active) values (?, ?, 1)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @process_data_ids;
  foreach my $pair (@$key_val_pairs) {
    $sth->bind_param(1, $pair->[0]);
    $sth->bind_param(2, $pair->[1]);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    push(@process_data_ids, $sth->{'mysql_insertid'});
  }
  return \@process_data_ids;
}

sub _add_branch_data {
  my ($self, $branch_ids, $key_val_pairs,) = @_;
  my $process_data_ids = $self->_add_process_data($key_val_pairs);
  my $sql = "INSERT INTO branch_data (branch_id, process_data_id) values (?, ?)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  foreach my $process_data_id (@$process_data_ids) {
    foreach my $branch_id (@$branches) {
      $sth->bind_param(1, $branch_id);
      $sth->bind_param(2, $process_data_id);
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    }
  }
}

sub _add_branches {
  my ($self, $num_branches, $branch_system_id) = @_;
  my $parent_branch_id = $self->{'_branch_id_hash'}->{$branch_system_id};
  my $sql = "INSERT INTO branch (parent_branch_id, creator_analysis_id, sibling_index, branch_system_id) VALUES (?, ?, ?, ?)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @new_branch_ids;
  foreach my $sibling_index (0..$num_branches-1) {
    $sth->bind_param(1, $parent_branch_id);
    $sth->bind_param(2, $self->analysis->dbID);
    $sth->bind_param(3, $sibling_index);
    $sth->bind_param(4, $branch_system_id);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    push(@new_branch_ids, $sth->{'mysql_insertid'});
  }
  return @new_branch_ids;
}









sub _branch_data {
  my ($self, %args) = @_;
  throw "no key" if ! $args{'key'};

  my $root_branch_id = $args{'branch_id'} || $self->param('branch_id');
  my @branch_ids = ($root_branch_id);
  if ($args{'descend'}) {
    push(@branch_ids, @{$self->_get_all_child_branch_ids($root_branch_id)});
  }
  if ($args{'ascend'}) {
    push(@branch_ids, @{$self->_get_all_parent_branch_ids($root_branch_id)});
  }
  my $hive_dbc = $self->data_dbc();

  if (defined $args{'value'}) {
    my $sql = "INSERT INTO process_data (data_key, data_value, is_active) values (?, ?, 1)";
    my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
    VALUE:
    foreach my $value (ref($args{'value'}) eq 'ARRAY' ? @{$args{'value'}} : ($args{'value'})) {
      next VALUE if ! defined $value;
      foreach my $branch_id (@branch_ids) {
        $sth->bind_param(1, $branch_id);
        $sth->bind_param(2, $args{'key'});
        $sth->bind_param(3, $value);
        $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      }
    }
  }
  else {
    my $sql = "SELECT branch_meta_data_id, meta_value FROM branch_meta_data WHERE branch_id=? AND meta_key=? AND is_active=1";
    my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
    my @return_data;
    foreach my $branch_id (@branch_ids) {
      $sth->bind_param(1, $branch_id);
      $sth->bind_param(2, $args{'key'});
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      my $rows = $sth->fetchall_hashref('branch_meta_data_id');
      push(@return_data, values %$rows);

      #my $rows = $sth->fetchall_arrayref;
      #push(@values, map {$_->[0]} @$rows);
    }
    return \@return_data;
  }
}

sub _get_all_child_branch_ids {
  my ($self, $parent_branch_id) = @_;
  $parent_branch_id ||= $self->param('branch_id');
  throw("no branch_id") if ! defined $parent_branch_id;
  my @child_branch_ids;
  my $hive_dbc = $self->data_dbc();
  my $sql = "SELECT child_branch_id FROM branch WHERE parent_branch_id=? ORDER BY child_branch_index";
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @unprocessed_branch_ids = ($parent_branch_id);
  while (@unprocessed_branch_ids) {
    my $branch_id = shift @unprocessed_branch_ids;
    $sth->bind_param(1, $branch_id);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    my $rows = $sth->fetchall_arrayref;
    push(@child_branch_ids, map {$_->[0]} @$rows);
    unshift(@unprocessed_branch_ids, map {$_->[0]} @$rows);
  }
  return \@child_branch_ids;
}

sub _get_all_parent_branch_ids {
  my ($self, $child_branch_id) = @_;
  $child_branch_id ||= $self->param('branch_id');
  throw("no branch_id") if ! defined $child_branch_id;
  my @parent_branch_ids;
  my $hive_dbc = $self->data_dbc();
  my $sql = "SELECT parent_branch_id FROM branch WHERE child_branch_id=?";
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my $current_branch_id = $child_branch_id;
  BRANCH:
  while (1) {
    $sth->bind_param(1, $current_branch_id);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    my $row = $sth->fetchrow_arrayref;
    last BRANCH if !defined $row || !defined $row->[0];
    $current_branch_id = $row->[0];
    push(@parent_branch_ids, $current_branch_id)
  }
  return \@parent_branch_ids;
}

sub _files_to_delete {
  my ($self, $arg) = @_;
  if ($arg) {
    $self->param('_files_to_delete', $self->_files_to_delete);
    my $files = ref($arg) eq 'ARRAY' ? $arg : [$arg];
    push(@{$self->param('_files_to_delete')}, @$files);
  }
  return $self->param('_files_to_delete') || [];
}

sub _data_to_make_inactive {
  my ($self, $arg) = @_;
  if ($arg) {
    $self->param('_data_to_make_inactive', $self->_data_to_make_inactive);
    my $dbIDs = ref($arg) eq 'ARRAY' ? $arg : [$arg];
    push(@{$self->param('_data_to_make_inactive')}, @$dbIDs);
  }
  return $self->param('_data_to_make_inactive') || [];
}

sub _inactivate_data {
  my ($self) = @_;
  my $sql = "UPDATE branch_meta_data SET is_active=0 WHERE branch_meta_data_id=?";
  my $dbIDs = $self->_data_to_make_inactive;
  return if ! scalar @$dbIDs;
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  foreach my $dbID (@$dbIDs) {
    $sth->bind_param(1, $dbID);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
  }
}


1;

