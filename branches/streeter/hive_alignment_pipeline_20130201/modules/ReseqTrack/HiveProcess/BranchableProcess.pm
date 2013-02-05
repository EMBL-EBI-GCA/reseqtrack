
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

  $self->_process_branch_params_in($branch_params_in);

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
    my @system_ids = keys %$branch_id_hash;

    my %output_data;
    while (my ($param_key, $value) = each %{$self->output_this_branch}) {
      my $db_hash = $branch_params_out->{$param_key};
      my $db_key = $db_hash ? $db_hash->{'key'} : $param_key;
      $output_data{$db_key} = $value;
    }
    $self->_add_branch_data(\@branch_ids, \%output_data);

    my $output_child_branches = $self->output_child_branches;
    my $branching_system = $self->param('branching_system');
    if (! defined $branching_system && scalar @system_ids == 1) {
      $branching_system = $system_ids[0];
    }
    my $num_child_branches = scalar @$output_child_branches;
    my $child_branch_ids = $self->_add_branches($num_child_branches, $branching_system);
    my @other_parents = map {$branch_id_hash->{$_}} grep {$_ != $branching_system} @system_ids;
    foreach my $i (0..$num_child_branches-1) {
      my @all_branch_ids = ($child_branch_ids->[$i], @other_parents);
      my %db_data;
      DATA_TYPE:
      while (my ($param_key, $value) = each %{$self->output_child_branches->[$i]}) {
        my $db_hash = $branch_params_out->{$param_key};
        my $db_key = $db_hash ? $db_hash->{'key'} : $param_key;
        $db_data{$db_key} = $value;
      }
      $self->_add_branch_data(\@all_branch_ids, \%db_data);
      $self->dataflow_output_id( {'branch_id' => \@all_branch_ids}, 2);
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
  throw("root_output_dir not defined") if !$root_output_dir;
  my $sub_dirs = $self->get_param_array('branch_subdir');
  my $output_dir = join('/', $root_output_dir, reverse @$sub_dirs);
  $self->output_dir($output_dir);
}
sub _process_labels {
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
    my $num_parent_ids = 0;
    my $query_branch_id = $branch_id;
    BRANCH:
    while (1) {
      $sth->bind_param(1, $query_branch_id);
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      my $row = $sth->fetchrow_arrayref;
      last BRANCH if !defined $row || !defined $row->[0];
      $query_branch_id = $row->[0];
      $parent_branch_ids{$system_id}{$query_branch_id} = $num_parent_ids;
      $num_parent_ids += 1;
    }
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
    my $ordered_branch_ids = __recursively_query_child_branch($branch_id, $sth);
    my $num_branch_ids = scalar @$ordered_branch_ids;
    foreach my $i (0..$num_branch_ids-1) {
      $child_branch_ids{$system_id}->{$ordered_branch_ids->[$i]} = $i;
    }
  }
  $self->{'_child_branch_id_hash'} = \%child_branch_ids;
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
  my $child_branch_id_hash = $self->{'_child_branch_id_hash'};
  my $parent_branch_id_hash = $self->{'_parent_branch_id_hash'};
  my $root_output_dir = $self->param('root_output_dir');

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
  PARAM:
  while (my ($param_key, $arg) = each %$branch_params_in) {
    my ($db_key, $ascend, $descend, $inactivate) = ref($arg) eq 'HASH'
                  ? @{$arg}{qw(key ascend descend inactivate)}
                  : ($arg, 0, 0, 0);
    my %data_value_hash;
    my %data_branch_hash;
    while (my ($system_id, $branch_id) = each %$branch_id_hash) {
      my @query_branch_ids = ($branch_id);
      push(@query_branch_ids, keys %{$parent_branch_id_hash->{$system_id}}) if $ascend;
      push(@query_branch_ids, keys %{$child_branch_id_hash->{$system_id}}) if $descend;
      foreach my $query_branch_id (@query_branch_ids) {
        $sth1->bind_param(1, $query_branch_id);
        $sth1->bind_param(2, $db_key);
        $sth1->execute() or die "could not execute $sql1: ".$sth1->errstr;
        VALUE:
        while (my $row1 = $sth1->fetchrow_arrayref()) {
          my ($data_id, $data_value, $num_systems) = @$row1;
          if ($num_systems == 1) {
            $data_value_hash{$data_id} = $data_value;
            $data_branch_hash{$data_id}{$system_id} = $query_branch_id;
            next VALUE;
          }
          next VALUE if defined $data_value_hash{$data_id};
          $sth2->bind_param(1, $data_id);
          $sth2->execute() or die "could not execute $sql2: ".$sth2->errstr;
          my $rows2 = $sth2->fetchall_arrayref();
          BRANCH:
          foreach my $row2 (@$rows2) {
            my ($other_branch_id, $other_system_id) = @$row2;
            next BRANCH if $branch_id_hash->{$other_system_id} == $other_branch_id;
            next BRANCH if $descend && defined $child_branch_id_hash->{$other_system_id}->{$other_branch_id};
            next BRANCH if $ascend && defined $parent_branch_id_hash->{$other_system_id}->{$other_branch_id};
            next VALUE;
          }
          $data_value_hash{$data_id} = $data_value;
          foreach my $row2 (@$rows2) {
            my ($other_branch_id, $other_system_id) = @$row2;
            $data_branch_hash{$data_id}{$other_system_id} = $other_branch_id;
          }
        }
      }
    }
    next PARAM if ! scalar keys %data_value_hash;
    #now to put the data in order.
    my @sorted_data_ids =
                  sort {__compare_branch_ids($a, $b, \%data_branch_hash, $child_branch_id_hash, $parent_branch_id_hash)}
                  keys %data_value_hash;
    my @sorted_values = map {$data_value_hash{$_}} @sorted_data_ids;
    $self->param($param_key, scalar @sorted_values == 1 ? $sorted_values[0] : \@sorted_values);

    if ($inactivate) {
      $self->_data_to_make_inactive(\@sorted_data_ids);
      $self->_files_to_delete([grep { /^$root_output_dir/ }  @sorted_values]);
    }
  }
}

sub __compare_branch_ids {
  my ($data_id_a, $data_id_b, $data_branch_hash, $child_branch_id_hash, $parent_branch_id_hash) = @_;
  my @systems_a = sort {$a <=> $b} keys %{$data_branch_hash->{$data_id_a}};
  my @systems_b = sort {$a <=> $b} keys %{$data_branch_hash->{$data_id_b}};
  my $ret = -(scalar @systems_a <=> scalar @systems_b);
  return $ret if $ret;
  foreach my $i (0..$#systems_a) {
    $ret = $systems_a[$i] <=> $systems_b[$i];
    return $ret if $ret;
  }


  # now look at parent systems
  SYSTEM:
  foreach my $system_id (@systems_a) {
    my $order_a = $parent_branch_id_hash->{$system_id}->{$data_branch_hash->{$data_id_a}{$system_id}};
    my $order_b = $parent_branch_id_hash->{$system_id}->{$data_branch_hash->{$data_id_b}{$system_id}};
    $ret = defined $order_a <=> defined $order_b;
    return $ret if $ret;
    next SYSTEM if ! defined $order_a;
    $ret = $order_a <=> $order_b;
    return $ret if $ret;
  }

  # now look at child systems
  SYSTEM:
  foreach my $system_id (@systems_a) {
    my $order_a = $child_branch_id_hash->{$system_id}->{$data_branch_hash->{$data_id_a}{$system_id}};
    my $order_b = $child_branch_id_hash->{$system_id}->{$data_branch_hash->{$data_id_b}{$system_id}};
    $ret = defined $order_a <=> defined $order_b;
    return $ret if $ret;
    next SYSTEM if ! defined $order_a;
    $ret = $order_a <=> $order_b;
    return $ret if $ret;
  }
  return 0;
}

sub get_param_array {
  my ($self, $key) = @_;
  my $param = $self->param($key);
  return $param if ref($param) eq 'ARRAY';
  return [] if !defined $param;
  return [$param];
}

sub output_child_branches {
  my ($self, %args) = @_;
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
  my ($self, $data_hash) = @_;
  my $sql = "INSERT INTO process_data (data_key, data_value, is_active) values (?, ?, 1)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @process_data_ids;
  while (my ($key, $vals) = each %$data_hash) {
    $vals = ref($vals) eq 'ARRAY' ? $vals : [$vals];
    foreach my $val (@$vals) {
      $sth->bind_param(1, $key);
      $sth->bind_param(2, $val);
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      push(@process_data_ids, $sth->{'mysql_insertid'});
    }
  }
  return \@process_data_ids;
}

sub _add_branch_data {
  my ($self, $branch_ids, $data_hash,) = @_;
  my $process_data_ids = $self->_add_process_data($data_hash);
  my $sql = "INSERT INTO branch_data (branch_id, process_data_id) values (?, ?)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  foreach my $process_data_id (@$process_data_ids) {
    foreach my $branch_id (@$branch_ids) {
      $sth->bind_param(1, $branch_id);
      $sth->bind_param(2, $process_data_id);
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    }
  }
}

sub _add_branches {
  my ($self, $num_branches, $branch_system_id) = @_;
  my $parent_branch_id = $self->{'_branch_id_hash'}->{$branch_system_id};
  my $sql = "INSERT INTO branch (parent_branch_id, sibling_index, branch_system_id) VALUES (?, ?, ?)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @new_branch_ids;
  foreach my $sibling_index (0..$num_branches-1) {
    $sth->bind_param(1, $parent_branch_id);
    $sth->bind_param(2, $sibling_index);
    $sth->bind_param(3, $branch_system_id);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    push(@new_branch_ids, $sth->{'mysql_insertid'});
  }
  return \@new_branch_ids;
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
  my $sql = "UPDATE process_data SET is_active=0 WHERE process_data_id=?";
  my $dbIDs = $self->_data_to_make_inactive;
  return if ! scalar @$dbIDs;
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  foreach my $dbID (@$dbIDs) {
    $sth->bind_param(1, $dbID);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
  }
}

sub output_dir {
  my ($self, $dir) = @_;
  if (defined $dir) {
    $self->{'output_dir'} = $dir;
  }
  return $self->{'output_dir'};
}

sub job_name {
  my ($self, $dir) = @_;
  if (defined $dir) {
    $self->{'job_name'} = $dir;
  }
  return $self->{'job_name'};
}


1;

