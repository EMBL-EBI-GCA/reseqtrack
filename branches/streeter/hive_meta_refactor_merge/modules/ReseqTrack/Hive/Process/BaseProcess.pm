=head1 NAME

ReseqTrack::Hive::Process::BaseProcess;

=head1 SYNOPSIS

    Extended Hive API for ReseqTrack pipelines

    Any ReseqTrack hive process should inherit from this base class.
    The fetch_input and write_output subroutines are defined here (requirements for any Hive Process)
    The child class should implement the run subroutine.

    Advantages of using this base class:

      Files can be stored in the reseqtrack_file table of the hive database.  A job input_id can then contain the file dbID instead of a long file string.
          (This behaviour can be turned off for a pipeline by setting 'use_reseqtrack_file_table' => 0 in the pipeline configuration file)

      Files can be deleted automatically after successful completion of a job. E.g. by adding the following in your analysis definition:
        -parameters => {
            delete_param => 'fastq'
        }
        or:
        -parameters => {
            delete_param => ['bam', 'bai']
        }
      Files will only be deleted if they are located somewhere under the root_output_dir

      Output files and directories are given sensible names.  This is managed by using the concept of a branch 'label'.  Each time the pipeline branches into factory jobs, each branch is given a new label string.
          (This behaviour can be turned off for a pipeline by setting 'use_label_management' => 0 in the pipeline configuration file)
          (Output directory can be overriden for a particular analysis by explicitly setting the parameter 'output_dir' => '/path/to/dir' in the pipeline configuration file)

      Management of job input_ids: When job A creates job B, the input_id for job B will contain:
          1. any parameter from the input_id of job A
          2. any parameter from the accu table for job A
          3. any extra data added by job A

=cut


package ReseqTrack::Hive::Process::BaseProcess;

use strict;
use ReseqTrack::Tools::Exception qw(throw);

use base ('Bio::EnsEMBL::Hive::Process');
use ReseqTrack::Tools::FileSystemUtils qw(delete_file);
use Bio::EnsEMBL::Hive::Utils qw(destringify);


=head2 output_param

  Arg [1]   : param_name, string
  Arg [2]   : arg, any Perl structure or object
  Function  : This data will be written to the dataflow_output_id when jobs are created later by the write_output subroutine
  Example   : $self->output_param('bam', '/path/to/new_file.bam');

=cut


sub output_param {
  my ($self, $param_name, $arg) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (@_ >=3 ) {
    $base_params->{'output_hash'}->{$param_name} = $arg;
  }
  if (@_ >=2) {
    return $base_params->{'output_hash'}->{$param_name};
  }
  return $base_params->{'output_hash'} // {};
}

=head2 prepare_factory_output_id

  Arg [1]   : hashref, any data unique to the child branch
  Function  : A factory module should create jobs by using this subroutine
              It should be called once for each child branch.  Each branch should have a different label and data hashref
              The jobs will be created later by the write_output subroutine
  Example   : $self->prepare_factory_output_id({run_id => 'ERR000001'});

=cut

sub prepare_factory_output_id {
  my ($self, $data_hash) = @_;
  throw('prepare_factory_output_id needs a data hash') if ref($data_hash) ne 'HASH';
  my $base_params = $self->param('_BaseProcess_params');
  push(@{$base_params->{'factory_outputs'}}, $data_hash);
}


=head2 output_dir

  Arg [1]   : (optional) dir
  Function  : Returns a sensible location to write all output for this job.
  Example   : my $output_dir = $self->output_dir
  Returntype: string

=cut


sub output_dir {
  my ($self, $dir) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (defined $dir) {
    $base_params->{'output_dir'} = $dir;
  }
  elsif($self->param_is_defined('root_output_dir')) {
    if ($self->param_is_defined('dir_label_params')) {
      $base_params->{'output_dir'} = join('/', $self->param('root_output_dir'),
              grep {length($_)}
              map {$self->param($_)}
              grep {$self->param_is_defined($_)}
              grep {length($_)}
              @{$self->param('dir_label_params')});
    }
    else {
      $base_params->{'output_dir'} = $self->param('root_output_dir');
    }
    $base_params->{'output_dir'} =~ s/\s+//g;
    $base_params->{'output_dir'} =~ s{//+}{/}g;
  }
  throw('cannot make a sensible output directory') if ! $base_params->{'output_dir'};
  return $base_params->{'output_dir'};
}

=head2 job_name

  Arg [1]   : (optional) job_name
  Function  : Returns a sensible job_name for passing to a RunProgram object. Unique within the pipeline.
  Example   : my $job_name = $self->job_name;
  Returntype: string

=cut

sub job_name {
  my ($self) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (!$base_params->{'job_name'}) {
    my $analysis_label = $self->param_is_defined('analysis_label') ?
              $self->param('analysis_label') : $self->analysis->logic_name;

    my $job_label;
    if ($self->param_is_defined('file_label')) {
      $job_label = $self->param('file_label');
    }
    if (!length($job_label) && $self->param_is_defined('dir_label_params')) {
      $job_label = ( grep {length($_)}
            map {$self->param($_)}
            grep {$self->param_is_defined($_)}
            grep {length($_)}
            @{$self->param('dir_label_params')})[-1];
    }
    if (!length($job_label)) {
      $job_label = $self->input_job->dbID;
    }
    $job_label =~ s/\s+//g;

    $base_params->{'job_name'} = join('.', $job_label, $analysis_label);
  }
  return $base_params->{'job_name'};
}


=head2 count_param_values

  Arg [1]   : param_name
  Function  : Counts how many values have been passed to the job for a particular parameter
            E.g. if $self->param('A') is [10,20,30] then $self->count_param_values('A') is 3
  Example   : my $num_bams = $self->count_param_values('bam');
  Returntype: Integer

=cut


sub count_param_values {
  my ($self, $param_name) = @_;
  return 0 if ! $self->param_is_defined($param_name);
  my $param_value = $self->param($param_name);
  my $flattened_values = [];
  __add_to_flat_array($flattened_values, $param_value);
  return scalar @$flattened_values;
}

=head2 file_param

  Arg [1]   : param_name
  Function  : Getter for job parameters.  Does a dbID to filename conversion if the 
            stored object is a file stored in the reseqtrack_file table (or a data structure
            containing files in the reseqtrack_file table)
  Example   : my $bam = $self->file_param('bam');
  Returntype: Any Perl structure or object

=cut

sub file_param {
  my ($self, $param_name) = @_;
  if (!$self->param('use_reseqtrack_file_table')) {
    return $self->param($param_name)
  }
  return undef if ! $self->param_is_defined($param_name);
  my $data_structure = $self->param($param_name);
  $data_structure = $self->_get_files($data_structure);
  return $data_structure;
}

=head2 param_to_flat_array

  Arg [1]   : param_name
  Function  : Getter for job parameters.  Converts complex data structures to arrays.
            E.g. a nested array [[[1,2],[4,5]],[7,8]] to [1,2,4,5,7,8]
  Returntype: Array ref
  Example   : my $bam_list = $self->file_param('bam');

=cut

sub param_to_flat_array {
  my ($self, $param_name) = @_;
  return [] if ! $self->param_is_defined($param_name);
  my $data_structure = $self->param($param_name);
  my $flattened_values = [];
  __add_to_flat_array($flattened_values, $data_structure);
  return $flattened_values;
}

=head2 file_param_to_flat_array

  Arg [1]   : param_name
  Function  : Getter for job parameters.  Combines the functions of param_to_flat_array and file_param.
  Returntype: Array ref
  Example   : my $bam_list = $self->file_param_to_flat_array('bam');

=cut


sub file_param_to_flat_array {
  my ($self, $param_name) = @_;
  my $flattened_values = $self->param_to_flat_array($param_name);
  if ($self->param('use_reseqtrack_file_table')) {
    $flattened_values = $self->_get_files($flattened_values);
  }
  return $flattened_values;
}

=head2 run_program

  Arg [1]   : A ReseqTrack::Tools::RunProgram object
  Arg [2]   : Additional arguments to pass to the RunProgram run subroutine
  Function  : A safe way to call the run subroutine of a RunProgram object in a hive analysis
  Exceptions: throws if execution of command fails
              Deliberately hangs if it notices a term_sig
  Example   : $self->run_program($run_samtools_object, 'merge');

=cut


sub run_program {
  my ($self, $run_program_object, @args) = @_;

  $self->dbc->disconnect_when_inactive(1);
  my $return = eval{$run_program_object->run(@args);};
  my $msg_thrown = $@;
  $self->dbc->disconnect_when_inactive(0);
  return $return if !$msg_thrown;
  print "term_sig is ".$run_program_object->term_sig . "\n";
  if ($run_program_object->term_sig) {
    while (1) {
      next; # Looks like sleep doesn't work when we have a term_sig
    }
  }
  die $msg_thrown;
}

=head2 fetch_input

  Description : All hive processes run three subroutines in order:
                    1. fetch_input
                    2. run
                    3. write_output
                This is the standard fetch_input subroutine for any ReseqTrack::Hive process
                Any ReseqTrack::Hive::Process class inheriting from this class only needs to define the run subroutine.

=cut

sub fetch_input {
  my ($self) = @_;

  $self->param_required('use_reseqtrack_file_table');

  $self->param('_BaseProcess_params', {});

  if ($self->param_is_defined('delete_param') && $self->param_is_defined('root_output_dir')) {
    my $root_output_dir = $self->param('root_output_dir');
    my %delete_params;
    my @delete_files;
    PARAM:
    foreach my $param_name (@{$self->param_to_flat_array('delete_param')}) {
      $delete_params{$param_name} = 1;
      push(@delete_files, @{$self->file_param_to_flat_array($param_name)});
    }
    $self->_files_to_delete([grep { /^$root_output_dir/ } @delete_files]);
    $self->_params_to_delete([keys %delete_params]);
  }

  my $flows_factory = $self->param_is_defined('flows_factory') ? $self->param('flows_factory') : 2;
  $self->_flows_factory($flows_factory);

  my $flows_non_factory = $self->param_is_defined('flows_non_factory') ? $self->param('flows_non_factory')
                        : (grep {$_ == 1} @{$self->_flows_factory}) ? undef : 1;
  $self->_flows_non_factory($flows_non_factory);

  my $flows_do_count = $self->param_is_defined('flows_do_count') ? $self->param('flows_do_count') : {};
  $self->_flows_do_count($flows_do_count);

  my $do_count_param = $self->param_is_defined('flows_do_count_param') ? $self->param('flows_do_count_param') : {};
  $self->_flows_do_count_param($do_count_param);

  throw("flows_do_count_param is not defined") if ! $do_count_param && scalar keys %$flows_do_count;


}

=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
}

=head2 write_output

  Description : All hive processes run three subroutines in order:
                    1. fetch_input
                    2. run
                    3. write_output
                This is the standard write_output subroutine for any ReseqTrack::Hive process
                Any ReseqTrack::Hive::Process class inheriting from this class only needs to define the run subroutine.

=cut


sub write_output {
  my ($self) = @_;

  $self->input_job->autoflow(0);

  my %base_output;
  my $use_file_table = $self->param('use_reseqtrack_file_table');

  my $delete_params = $self->_params_to_delete;
  foreach my $param (@$delete_params) {
    $base_output{$param} = undef;
  }

  while (my ($key, $val) = each %{$self->output_param}) {
    if ($use_file_table) {
      $val = $self->_make_file_ids($val);
    }
    $self->param($key, $val);
    $base_output{$key} = $val;
  }

  my $flows_factory = $self->_flows_factory();
  my $flows_non_factory = $self->_flows_non_factory();
  my $flows_do_count_hash = $self->_flows_do_count();
  my $flows_do_count_param = $self->_flows_do_count_param();

  my $num_files_base = $flows_do_count_param ? $self->count_param_values($flows_do_count_param) : 0;
  my %flows_pass_count_base;
  FLOW:
  foreach my $flow (@$flows_non_factory, @$flows_factory) {
    if ($flows_do_count_hash->{$flow}) {
      my ($require_num, $modifier) = @{$flows_do_count_hash->{$flow}};
      next FLOW if !$modifier && $require_num != $num_files_base;
      next FLOW if $modifier eq '+' && $require_num > $num_files_base;
      next FLOW if $modifier eq '-' && $require_num < $num_files_base;
    }
    $flows_pass_count_base{$flow} = 1;
  }


  my $factory_outputs = $self->param('_BaseProcess_params')->{'factory_outputs'};
  if (defined $factory_outputs && @$factory_outputs && @$flows_factory) {
    my %flows_outputs;
    my @flows_for_count = grep {defined $flows_do_count_hash->{$_}} @$flows_factory;
    my @flows_no_count = grep {! defined $flows_do_count_hash->{$_}} @$flows_factory;
    OUTPUT:
    foreach my $extra_data_hash (@$factory_outputs) {
      my @flows_for_output = @flows_no_count;
      if (@flows_for_count) {
        if (exists $extra_data_hash->{$flows_do_count_param}) {
          my $flattened_values = [];
          __add_to_flat_array($flattened_values, $extra_data_hash->{$flows_do_count_param});
          my $num_files = scalar @$flattened_values;
          FLOW:
          foreach my $flow (@flows_for_count) {
            my ($require_num, $modifier) = @{$flows_do_count_hash->{$flow}};
            next FLOW if !$modifier && $require_num != $num_files;
            next FLOW if $modifier eq '+' && $require_num > $num_files;
            next FLOW if $modifier eq '-' && $require_num < $num_files;
            push(@flows_for_output, $flow);
          }
        }
        else {
          push(@flows_for_output, grep {$flows_pass_count_base{$_}} @flows_for_count);
        }
      }
      next OUTPUT if !@flows_for_output;

      my $new_factory_hash = {%base_output};
      while (my ($param_name, $param_value) = each %$extra_data_hash) {
        if ($use_file_table) {
          $param_value = $self->_make_file_ids($param_value);
        }
        $new_factory_hash->{$param_name} = $param_value;
      }
      foreach my $flow (@flows_for_output) {
        push(@{$flows_outputs{$flow}}, $new_factory_hash);
      }
    }
    while (my ($flow, $output_hashes) = each %flows_outputs) {
      $self->dataflow_output_id($output_hashes, $flow);
    }
  }

  foreach my $flow (sort {$b <=> $a} grep {$flows_pass_count_base{$_}} @$flows_non_factory) {
    $self->dataflow_output_id(\%base_output, $flow);
  }

  foreach my $file (@{$self->_files_to_delete}) {
    delete_file($file);
  }
}

# ########## private subroutines that should not be called by child class ###############
sub _flows_non_factory {
  my ($self, $arg) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (@_ >= 2) {
    if (!defined $arg) {
      $base_params->{'flows_non_factory'} = [];
    }
    elsif (ref($arg) eq 'ARRAY') {
      $base_params->{'flows_non_factory'} = $arg;
    }
    elsif (ref($arg) eq 'HASH') {
      my @flows;
      while (my ($flow, $is_true) = each %$arg) {
        push(@flows, $flow) if $is_true;
      }
      $base_params->{'flows_non_factory'} = \@flows;
    }
    else {
      $base_params->{'flows_non_factory'} = [$arg];
    }
  }
  return $base_params->{'flows_non_factory'} // [];
}

sub _flows_factory {
  my ($self, $arg) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (@_ >= 2) {
    if (!defined $arg) {
      $base_params->{'flows_factory'} = [];
    }
    elsif (ref($arg) eq 'ARRAY') {
      $base_params->{'flows_factory'} = $arg;
    }
    elsif (ref($arg) eq 'HASH') {
      my @flows;
      while (my ($flow, $is_true) = each %$arg) {
        push(@flows, $flow) if $is_true;
      }
      $base_params->{'flows_factory'} = \@flows;
    }
    else {
      $base_params->{'flows_factory'} = [$arg];
    }
  }
  return $base_params->{'flows_factory'} // [];
}

sub _flows_do_count {
  my ($self, $hash) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (@_ >= 2) {
    throw('not a hash') if ref($hash) ne 'HASH';
    while (my ($flow, $condition) = each %$hash) {
      my ($require_num, $modifier) = $condition =~ /(\d+)([+-]?)/;
      throw("did not recognise condition for $flow") if !defined $require_num;
      $base_params->{'flows_do_count'}->{$flow} = [$require_num, $modifier];
    }
  }
  return $base_params->{'flows_do_count'} // {};
}

sub _flows_do_count_param {
  my ($self, $param) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (defined $param) {
    $base_params->{'flows_do_count_param'} = $param;
  }
  return $base_params->{'flows_do_count_param'};
}



sub _get_files {
  my ($self, $data_structure) = @_;
  my $sql = "SELECT name FROM reseqtrack_file WHERE file_id=?";
  #my $hive_dbc = $self->data_dbc();
  my $hive_dbc = $self->dbc();
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  return $self->__structure_to_file_paths($data_structure, $sth);
}

sub _make_file_ids {
  my ($self, $data_structure) = @_;
  #my $hive_dbc = $self->data_dbc();
  my $hive_dbc = $self->dbc();
  my $sql_insert = "INSERT INTO reseqtrack_file (name) values (?)";
  my $sth_insert = $hive_dbc->prepare($sql_insert) or die "could not prepare $sql_insert: ".$hive_dbc->errstr;

  return __structure_to_file_ids($data_structure, $sth_insert);
}


sub __structure_to_file_ids {
  my ($data_structure, $sth_insert) = @_;
  if (ref($data_structure) eq 'ARRAY') {
    foreach my $i (0..$#{$data_structure}) {
      $data_structure->[$i] = __structure_to_file_ids($data_structure->[$i], $sth_insert);
    }
    return $data_structure;
  }
  elsif (ref($data_structure) eq 'HASH') {
    my %new_hash;
    while (my ($key, $value) = each %$data_structure) {
      my $new_key = __structure_to_file_ids($key, $sth_insert);
      my $new_value = __structure_to_file_ids($value, $sth_insert);
      $new_hash{$new_key} = $new_value;
    }
    return \%new_hash;
  }
  elsif ($data_structure =~ m{/\S}) {
    $sth_insert->bind_param(1, $data_structure);
    $sth_insert->execute() or die 'could not execute '.$sth_insert->statement .': '.$sth_insert->errstr;
    my $file_id = $sth_insert->{'mysql_insertid'};
    return 'F'.$file_id.'F';
  }
  else {
    return $data_structure;
  }
}

sub __structure_to_file_paths {
  my ($self, $data_structure, $sth) = @_;
  if (ref($data_structure) eq 'ARRAY') {
    foreach my $i (0..$#{$data_structure}) {
      $data_structure->[$i] = $self->__structure_to_file_paths($data_structure->[$i], $sth);
    }
    return $data_structure;
  }
  elsif (ref($data_structure) eq 'HASH') {
    my %new_hash;
    while (my ($key, $value) = each %$data_structure) {
      my $new_key = $self->__structure_to_file_paths($key, $sth);
      my $new_value = $self->__structure_to_file_paths($value, $sth);
      $new_hash{$new_key} = $new_value;
    }
    return \%new_hash;
  }
  elsif ($data_structure =~ /^F(\d+)F$/ ) {
    my $file_id = $1;
    if (my $file = $self->_file_cache($file_id)) {
      return $file;
    }
    else {
      $sth->bind_param(1, $file_id);
      $sth->execute() or die 'could not execute '.$sth->statement .': '.$sth->errstr;
      my $file = $sth->fetchall_arrayref->[0]->[0];
      return $data_structure if !$file;
      $self->_file_cache($file_id, $file);
      return $file;
    }
  }
  else {
    return $data_structure;
  }
}

sub _params_to_delete {
  my ($self, $arg) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (@_ >=2) {
    $base_params->{'params_to_delete'} = $arg;
  }
  return $base_params->{'params_to_delete'} // [];
}

sub _files_to_delete {
  my ($self, $arg) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (@_ >=2) {
    my $files = ref($arg) eq 'ARRAY' ? $arg : [$arg];
    push(@{$base_params->{'files_to_delete'}}, @$files);
  }
  return $base_params->{'files_to_delete'} // [];
}

sub _file_cache {
  my ($self, $dbID, $name) = @_;
  my $base_params = $self->param('_BaseProcess_params');
  if (defined $name) {
    $base_params->{'file_cache'}->{$dbID} = $name;
  }
  return $base_params->{'file_cache'}->{$dbID};
}

sub __add_to_flat_array {
  my ($flat_array, $value) = @_;
  if (ref($value) eq 'ARRAY') {
    foreach my $sub_value (@$value) {
      __add_to_flat_array($flat_array, $sub_value);
    }
  }
  elsif (defined $value) {
    push(@$flat_array, $value);
  }
}

1;

