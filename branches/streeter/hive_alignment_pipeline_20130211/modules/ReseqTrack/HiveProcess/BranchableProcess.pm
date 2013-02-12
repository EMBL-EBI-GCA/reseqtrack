
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

    my $branch_id = $self->_branch_id;

    my %output_data;
    while (my ($param_key, $value) = each %{$self->output_this_branch}) {
      my $db_hash = $branch_params_out->{$param_key};
      my $db_key = $db_hash ? $db_hash->{'key'} : $param_key;
      $output_data{$db_key} = $value;
    }
    $self->_add_branch_data($branch_id, \%output_data);

    my $output_child_branches = $self->output_child_branches;
    my $num_child_branches = scalar @$output_child_branches;
    my $child_branch_ids = $self->_add_branches($num_child_branches);
    foreach my $i (0..$num_child_branches-1) {
      my %db_data;
      DATA_TYPE:
      while (my ($param_key, $value) = each %{$self->output_child_branches->[$i]}) {
        my $db_hash = $branch_params_out->{$param_key};
        my $db_key = $db_hash ? $db_hash->{'key'} : $param_key;
        $db_data{$db_key} = $value;
      }
      $self->_add_branch_data($child_branch_ids->[$i], \%db_data);
      $self->dataflow_output_id( {'branch_id' => $child_branch_ids->[$i]}, 2);
    }
    
    foreach my $i (@{$self->flows_this_branch}) {
      $self->dataflow_output_id( {'branch_id' => $branch_id}, $i);
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
  $self->_branch_id($self->param('branch_id'));
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
  my $branch_id = $self->_branch_id;
  if (!defined $branch_id) {
    $self->{'_parent_branch_id_array'} = [];
    return;
  }
  my @parent_branch_ids;
  my $sql = "SELECT parent_branch_id FROM branch WHERE branch_id=?";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my $query_branch_id = $branch_id;
  BRANCH:
  while (1) {
    $sth->bind_param(1, $query_branch_id);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    my $row = $sth->fetchrow_arrayref;
    last BRANCH if !defined $row || !defined $row->[0];
    $query_branch_id = $row->[0];
    push(@parent_branch_ids, $query_branch_id);
  }
  $self->{'_parent_branch_id_hash'} = \@parent_branch_ids;
}

sub _process_child_branch_ids {
  my ($self) = @_;
  my $branch_id = $self->_branch_id;
  if (!defined $branch_id) {
    $self->{'_child_branch_id_array'} = [];
    return;
  }
  my $sql = "SELECT branch_id FROM branch WHERE parent_branch_id=? ORDER BY sibling_index";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my $ordered_branch_ids = __recursively_query_child_branch($branch_id, $sth);
  $self->{'_child_branch_id_array'} = $ordered_branch_ids;
}

sub __recursively_query_child_branch {
  my ($branch_id, $sth) = @_;
  my @child_branch_ids;
  $sth->bind_param(1, $branch_id);
  $sth->execute() or die "could not execute: ".$sth->errstr;
  my $rows = $sth->fetchall_arrayref();
  foreach my $child_branch_id (map {$_->[0]} @$rows) {
    push(@child_branch_ids, $child_branch_id, @{__recursively_query_child_branch($child_branch_id, $sth)});
  }
  return \@child_branch_ids;
}

sub _process_branch_params_in {
  my ($self, $branch_params_in) = @_;
  my $branch_id = $self->_branch_id;
  return if ! defined $branch_id;
  my $child_branch_id_array = $self->{'_child_branch_id_array'};
  my $parent_branch_id_array = $self->{'_parent_branch_id_hash'};
  my $root_output_dir = $self->param('root_output_dir');

  my $sql = "SELECT branch_data_id, data_value FROM branch_data WHERE branch_id=?  AND data_key = ?  AND is_active=1";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  PARAM:
  while (my ($param_key, $arg) = each %$branch_params_in) {
    my ($db_key, $ascend, $descend, $inactivate) = ref($arg) eq 'HASH'
                  ? @{$arg}{qw(key ascend descend inactivate)}
                  : ($arg, 0, 0, 0);
    my @other_branch_arrays;
    push(@other_branch_arrays, $child_branch_id_array) if $descend;
    push(@other_branch_arrays, $parent_branch_id_array) if $ascend;
    my @sorted_values;
    my @data_ids;
    foreach my $query_branch_id ($branch_id, map {@$_} @other_branch_arrays) {
      $sth->bind_param(1, $query_branch_id);
      $sth->bind_param(2, $db_key);
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      my $rows = $sth->fetchall_arrayref;
      push(@data_ids, map {$_->[0]} @$rows);
      push(@sorted_values, map {$_->[1]} @$rows);
    }
    $self->param($param_key, scalar @sorted_values == 1 ? $sorted_values[0] : \@sorted_values);
    if ($inactivate) {
      $self->_data_to_make_inactive(\@data_ids);
      $self->_files_to_delete([grep { /^$root_output_dir/ }  @sorted_values]);
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
    $self->{'_flows_this_branch'} = $flows;
  }
  return $self->{'_flows_this_branch'} || [1];
}

sub _add_branch_data {
  my ($self, $branch_id, $data_hash) = @_;
  my $sql = "INSERT INTO branch_data (branch_id, data_key, data_value, is_active) values (?, ?, ?, 1)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @branch_data_ids;
  while (my ($key, $vals) = each %$data_hash) {
    $vals = ref($vals) eq 'ARRAY' ? $vals : [$vals];
    foreach my $val (@$vals) {
      $sth->bind_param(1, $branch_id);
      $sth->bind_param(2, $key);
      $sth->bind_param(3, $val);
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      push(@branch_data_ids, $sth->{'mysql_insertid'});
    }
  }
  return \@branch_data_ids;
}


sub _add_branches {
  my ($self, $num_branches) = @_;
  my $parent_branch_id = $self->_branch_id;
  my $sql = "INSERT INTO branch (parent_branch_id, sibling_index) VALUES (?, ?)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @new_branch_ids;
  foreach my $sibling_index (0..$num_branches-1) {
    $sth->bind_param(1, $parent_branch_id);
    $sth->bind_param(2, $sibling_index);
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
  my $sql = "UPDATE branch_data SET is_active=0 WHERE branch_data_id=?";
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

sub _branch_id {
  my ($self, $branch_id) = @_;
  if (defined $branch_id) {
    $self->{'_branch_id'} = $branch_id;
  }
  return $self->{'_branch_id'};
}


1;

