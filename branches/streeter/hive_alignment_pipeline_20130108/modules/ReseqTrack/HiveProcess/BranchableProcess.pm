
package ReseqTrack::HiveProcess::BranchableProcess;

use strict;
use ReseqTrack::Tools::Exception qw(throw);

use base ('Bio::EnsEMBL::Hive::Process');

=head2 fetch_input

    Description : Implements fetch_input() interface method of Bio::EnsEMBL::Hive::Process that is used to read in parameters and load data.
                  Here we have nothing to fetch.

=cut

sub fetch_input {
  my ($self) = @_;
  throw("do not have a branch id") if !defined $self->param('branch_id');
  my %branch_params = %{$self->param('universal_branch_parameters')};
  map {$branch_params{$_} = $self->param('branch_parameters')->{$_}} keys %{$self->param('branch_parameters')};
  while (my ($param_key, $arg) = each %branch_params) {
    my $hashref = ref($arg) eq 'HASH' ? $arg : {key => $arg};
    my @param_vals;
    my $db_keys = ref($hashref->{'key'}) eq 'ARRAY' ? $hashref->{'key'} : [$hashref->{'key'}];
    foreach my $db_key (@$db_keys) {
      push(@param_vals, @{$self->_branch_meta_data(key => $db_key, descend => $hashref->{'descend'}, ascend => $hashref->{'ascend'})});
    }
    if ($hashref->{'arrayref'}) {
      $self->param($param_key, \@param_vals);
    }
    else {
      throw("Found more than one value for $param_key in branch ".$self->param('branch_id')) if @param_vals > 1;
      $self->param($param_key, $param_vals[0]);
    }
  }

  if (!$self->param('output_dir')) {
    $self->param('output_dir', $self->param('root_output_dir'));
  }
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

    my $output_child_branches = $self->output_child_branches;

    my $num_child_branches = scalar @$output_child_branches;
    my $child_branch_ids = $self->_make_child_branches($num_child_branches);
    foreach my $i (0..$num_child_branches-1) {
      while (my ($key, $value) = each %{$output_child_branches->[$i]}) {
        $self->_branch_meta_data('branch_id' => $child_branch_ids->[$i],
                                'key' => $key, 'value' => $value);
      }
      $self->dataflow_output_id( {'branch_id' => $child_branch_ids->[$i]}, 2);
    }

    while (my ($key, $value) = each %{$self->output_this_branch}) {
      $self->_branch_meta_data('key' => $key, 'value' => $value);
    }
    $self->dataflow_output_id( {'branch_id' => $self->param('branch_id')}, 1);
}

sub output_child_branches {
  my ($self, %args) = @_;
  if (scalar keys %args) {
    $self->param('output_child_branches', $self->output_child_branches);
    push(@{$self->param('output_child_branches')}, \%args);
  }
  return $self->param('output_child_branches') || [];
}

sub output_this_branch {
  my ($self, %args) = @_;
  if (scalar keys %args) {
    $self->param('output_this_branch', \%args);
  }
  return $self->param('output_this_branch') || {};
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
    my $sql = "insert into branch_meta_data (branch_id, meta_key, meta_value) values (?, ?, ?)";
    my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
    foreach my $value (ref($args{'value'}) eq 'ARRAY' ? @{$args{'value'}} : ($args{'value'})) {
      foreach my $branch_id (@branch_ids) {
        $sth->bind_param(1, $branch_id);
        $sth->bind_param(2, $args{'key'});
        $sth->bind_param(3, $value);
        $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      }
    }
  }
  else {
    my $sql = "SELECT meta_value FROM branch_meta_data WHERE branch_id=? AND meta_key=?";
    my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
    my @values;
    foreach my $branch_id (@branch_ids) {
      $sth->bind_param(1, $branch_id);
      $sth->bind_param(2, $args{'key'});
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      my $rows = $sth->fetchall_arrayref;
      push(@values, map {$_->[0]} @$rows);
    }
    return \@values;
  }
}

sub _get_all_child_branch_ids {
  my ($self, $parent_branch_id) = @_;
  $parent_branch_id ||= $self->param('branch_id');
  throw("no branch_id") if ! defined $parent_branch_id;
  my @child_branch_ids;
  my $hive_dbc = $self->data_dbc();
  my $sql = "SELECT child_branch_id FROM branch WHERE parent_branch_id=?";
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @unprocessed_branch_ids = ($parent_branch_id);
  while (@unprocessed_branch_ids) {
    my $branch_id = shift @unprocessed_branch_ids;
    $sth->bind_param(1, $branch_id);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    my $rows = $sth->fetchall_arrayref;
    push(@child_branch_ids, map {$_->[0]} @$rows);
    push(@unprocessed_branch_ids, map {$_->[0]} @$rows);
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


1;

