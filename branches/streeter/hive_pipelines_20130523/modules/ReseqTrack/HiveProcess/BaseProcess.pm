
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

  my $root_output_dir = $self->param_required('root_output_dir');

  my $analysis_label = $self->param_is_defined('analysis_label') ?
            $self->param('analysis_label') : $self->analysis->logic_name;
  if ($self->param_is_defined('labels') {
    my $labels = $self->param('labels');
    throw("labels is defined but empty") if !@$labels;
    my $output_dir = join('/', $root_output_dir, @$labels};
    $output_dir =~ s{//}{/}g;
    $self->output_dir($output_dir);
    my $job_name = join('.', $labels->[-1], $analysis_label);
    $self->job_name($job_name);
  }
  else {
    $self->output_dir($root_output_dir);
    my $job_name = join('.', $self->input_job->dbID, $analysis_label);
    $self->job_name($job_name);
  }

  if ($self->param_is_defined('delete_param')) {
    my @file_ids;
    PARAM:
    foreach my $param_name (@{$self->get_param_values('delete_param')}) {
      push(@file_ids, @{$self->get_param_values($param_name)});
    }
    my $files = $self->get_files(\@file_ids);
    $self->_files_to_delete([grep { /^$root_output_dir/ } @$files);
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

  my $input_id_hash = destringify($self->input_id);
  my $params_to_add = $self->get_param_values('add_params_to_output');

  my %output_hash;
  foreach my $param_name (keys %$input_id_hash, @$params_to_add) {
    $output_hash{$param_name} = $self->param($param_name);
  }

  my $input_labels = $self->get_param_values('labels');
  foreach my $child (@{$self->{'_child_outputs'}}) {
    my ($child_label, $child_data_hash) = @$child;
    my @output_labels = (@$input_labels, $child_label);
    my %child_output_hash = %output_hash;
    while (my ($param_name, $param_value) = each %$child_data_hash) {
      $child_output_hash($param_name) = $param_value;
    }
    $child_output_hash('labels') = \@output_labels;
    $self->dataflow_output_id(\%child_output_hash, 2);
  }

  $self->dataflow_output_id(\%output_hash, 1);

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

sub get_files {
  my ($self, $file_ids) = @_;
  my $sql = "SELECT name FROM reseqtrack_file WHERE file_id=?";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @files;
  foreach my $file_id (@$file_ids) {
    if (my $file = $self->_file_cache($file_id)) {
      push(@files, $file);
    }
    else {
      $sth->bind_param(1, $file_id);
      $sth->execute() or die "could not execute $sql: ".$sth->errstr;
      my $rows = $sth->fetchall_arrayref;
      throw("nothing in database for file id $file_id") if !@$rows;
      $file = $rows->[0]->[0];
      $self->_file_cache($file_id, $file);
      push(@files, $file);
    }
  }
  return \@files;
}

sub make_file_ids {
  my ($self, $files) = @_;
  my $sql = "INSERT INTO reseqtrack_file (name) values (?)";
  my $hive_dbc = $self->data_dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  my @file_ids;
  foreach my $file (@$files) {
    $sth->bind_param(1, $file);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    push(@file_ids, $sth->{'mysql_insertid'});
  }
  return \@file_ids;
}

sub add_to_dataflow {
  my ($self, @param_names);
  my $added_params = $self->get_param_values('add_params_to_output');
  push(@$added_params, @param_names);
  $self->param('add_params_to_output', $added_params);
}

sub prepare_child_output_id {
  my ($self, $label, $data_hash) = @_;
  push(@{$self->{'_child_outputs'}}, [$label, $data_hash || {}]);
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


sub _file_cache {
  my ($self, $dbID, $name) = @_;
  if (defined $name) {
    $self->{'_file_cache'}->{$dbID} = $name;
  }
  return $self->{'_file_cache'}->{$dbID};
}

sub get_param_values {
  my ($self, $param_name) = @_;
  return [] if ! $self->param_is_defined($param_name);
  my $param_value = $self->param($param_name);
  my @flattened_values;
  __add_to_flat_array(\@flattened_values, $param_value);
  return \@flattened_values;
}
sub __add_to_flat_array {
  my ($flat_array, $value) = @_;
  if (ref($value) eq 'ARRAY') {
    foreach my $sub_value (@$value) {
      __add_to_flat_array($flat_array, $sub_value);
    }
  }
  else {
    push(@$flat_array, $value);
  }
}


1;

