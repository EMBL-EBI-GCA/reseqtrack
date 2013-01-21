
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
  throw("do not have a branch id") if !defined $self->param('branch_id');
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
        foreach my $file ( map {$_->{'meta_value'}} grep {! $_->{'never_delete'}} @$branch_data) {
          next if $file !~ /^$root_output_dir/;
          $self->_files_to_delete($file);
        }
      }
    }
    if (scalar @param_vals) {
      $self->param($param_key, scalar @param_vals == 1 ? $param_vals[0] : \@param_vals);
    }
  }

  if (!$self->param('output_dir')) {
    $self->param('output_dir', $root_output_dir);
  }

  my $analysis_name = $self->analysis->logic_name;
  my $branch_label = $self->param('branch_label');
  my $job_name = $branch_label ? "$branch_label.$analysis_name" : $analysis_name;
  $self->param('job_name', $job_name);
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

    DATA_TYPE:
    while (my ($param_key, $value) = each %{$self->output_this_branch}) {
      my $db_hash = $branch_params_out->{$param_key};
      if ($db_hash) {
        $self->_branch_meta_data(key => $db_hash->{'key'}, value => $value, never_delete => $db_hash->{'never_delete'});
      }
      else {
        $self->_branch_meta_data(key => $param_key, value => $value);
      }
    }

    my $output_child_branches = $self->output_child_branches;
    my $num_child_branches = scalar @$output_child_branches;
    my $child_branch_ids = $self->_make_child_branches($num_child_branches);
    foreach my $i (0..$num_child_branches-1) {
      DATA_TYPE:
      while (my ($param_key, $value) = each %{$self->output_child_branches->[$i]}) {
        my $db_hash = $branch_params_out->{$param_key};
        if ($db_hash) {
          $self->_branch_meta_data('branch_id' => $child_branch_ids->[$i],
                                'key' => $db_hash->{'key'}, 'value' => $value, never_delete => $db_hash->{'never_delete'});
        }
        else {
          $self->_branch_meta_data('branch_id' => $child_branch_ids->[$i], 'key' => $param_key, 'value' => $value);
        }
      }
      $self->dataflow_output_id( {'branch_id' => $child_branch_ids->[$i]}, 2);
    }

    foreach my $i (@{$self->flows_this_branch}) {
      $self->dataflow_output_id( {'branch_id' => $self->param('branch_id')}, $i);
    }

    foreach my $file (@{$self->_files_to_delete}) {
      delete_file($file);
    }
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
    $self->param('_flows_this_branch', $arg);
  }
  return $self->param('_flows_this_branch') || [1];
}

sub _make_child_branches {
  my ($self, $num_branches) = @_;

  my $sql = "insert into branch (parent_branch_id, child_branch_index) values (?, ?)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @child_branch_ids;
  foreach my $child_branch_index (0..$num_branches-1) {
    $sth->bind_param(1, $self->param('branch_id'));
    $sth->bind_param(2, $child_branch_index);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    push(@child_branch_ids, $sth->{'mysql_insertid'});
  }
  return \@child_branch_ids;
}

sub _branch_meta_data {
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
    my $sql = "insert into branch_meta_data (branch_id, meta_key, meta_value, is_active, never_delete values (?, ?, ?, 1, ?)";
    my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
    foreach my $value (ref($args{'value'}) eq 'ARRAY' ? @{$args{'value'}} : ($args{'value'})) {
      foreach my $branch_id (@branch_ids) {
        $sth->bind_param(1, $branch_id);
        $sth->bind_param(2, $args{'key'});
        $sth->bind_param(3, $value);
        $sth->bind_param(4, $args{'never_delete'} || 0);
        $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      }
    }
  }
  else {
    my $sql = "SELECT branch_meta_data_id, meta_value, never_delete FROM branch_meta_data WHERE branch_id=? AND meta_key=? AND is_active=1";
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

