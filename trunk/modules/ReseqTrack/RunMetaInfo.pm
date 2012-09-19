package ReseqTrack::RunMetaInfo;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::HasHistory;

@ISA = qw(ReseqTrack::HasHistory);


sub new{
  my ($class, @args) = @_;
  my  $self = $class->SUPER::new(@args);
  my ($run_id, $study_id, $study_name, $center_name, $submission_id,
      $submission_date, $sample_id, $sample_name, $population,
      $experiment_id, $instrument_platform, $instrument_model,
      $library_name, $run_name, $run_block_name, $paired_length,
      $status, $library_layout, $archive_base_count, $archive_read_count, $library_strategy) = 
      rearrange(['RUN_ID', 'STUDY_ID', 'STUDY_NAME', 'CENTER_NAME',
                 'SUBMISSION_ID', 'SUBMISSION_DATE', 'SAMPLE_ID', 'SAMPLE_NAME',
                 'POPULATION', 'EXPERIMENT_ID', 'INSTRUMENT_PLATFORM', 
                 'INSTRUMENT_MODEL', 'LIBRARY_NAME', 'RUN_NAME', 'RUN_BLOCK_NAME',
                 'PAIRED_LENGTH', 'STATUS', 'LIBRARY_LAYOUT', 'ARCHIVE_BASE_COUNT',
                 'ARCHIVE_READ_COUNT', 'LIBRARY_STRATEGY'], @args);
  
  #ERROR CHECKING
  throw("ReseqTrack::RunMetaInfo must have a run_id") unless($run_id);
  ##############
  $self->run_id($run_id);
  $self->study_id($study_id);
  $self->study_name($study_name);
  $self->center_name($center_name);
  $self->submission_id($submission_id);
  $self->submission_date($submission_date);
  $self->sample_id($sample_id);
  $self->sample_name($sample_name);
  $self->population($population);
  $self->experiment_id($experiment_id);
  $self->instrument_platform($instrument_platform);
  $self->instrument_model($instrument_model);
  $self->library_name($library_name);
  $self->run_name($run_name);
  $self->run_block_name($run_block_name);
  $self->paired_length($paired_length) if(defined($paired_length));
  $self->paired_length(0) unless(defined($self->paired_length));
  $self->library_layout($library_layout);
  $self->archive_base_count($archive_base_count) if(defined($archive_base_count));
  $self->archive_read_count($archive_read_count) if(defined($archive_read_count));
  $self->archive_base_count(0) unless(defined($self->archive_base_count));
  $self->archive_read_count(0) unless(defined($self->archive_read_count));
  $self->status($status);
  $self->library_strategy($library_strategy);
  $self->fix_date;
  return $self;
}

sub name{
  my ($self, $arg) = @_;
  if($arg){
    $self->run_id($arg);
  }
  return $self->run_id;
}

sub run_id{
  my ($self, $run_id) = @_;
  if($run_id){
    $self->{run_id} = $run_id;
  }
  return $self->{run_id};
}

sub study_id{
  my ($self, $study_id) = @_;
  if($study_id){
    $self->{study_id} = $study_id;
  }
  return $self->{study_id};
}

sub study_name{
  my ($self, $study_name) = @_;
  if($study_name){
    $self->{study_name} = $study_name;
  }
  return $self->{study_name};
}

sub center_name{
  my ($self, $center_name) = @_;
  if($center_name){
    $self->{center_name} = $center_name;
  }
  return $self->{center_name};
}

sub submission_id{
  my ($self, $submission_id) = @_;
  if($submission_id){
    $self->{submission_id} = $submission_id;
  }
  return $self->{submission_id};
}

sub submission_date{
  my ($self, $submission_date) = @_;
  if($submission_date){
    $self->{submission_date} = $submission_date;
  }
  return $self->{submission_date};
}

sub sample_id{
  my ($self, $sample_id) = @_;
  if($sample_id){
    $self->{sample_id} = $sample_id;
  }
  return $self->{sample_id};
}

sub sample_name{
  my ($self, $sample_name) = @_;
  if($sample_name){
    $self->{sample_name} = $sample_name;
  }
  return $self->{sample_name};
}

sub population{
  my ($self, $population) = @_;
  if($population){
    $self->{population} = $population;
  }
  return $self->{population};
}

sub experiment_id{
  my ($self, $experiment_id) = @_;
  if($experiment_id){
    $self->{experiment_id} = $experiment_id;
  }
  return $self->{experiment_id};
}

sub instrument_platform{
  my ($self, $instrument_platform) = @_;
  if($instrument_platform){
    $self->{instrument_platform} = $instrument_platform;
  }
  return $self->{instrument_platform};
}

sub instrument_model{
  my ($self, $instrument_model) = @_;
  if($instrument_model){
    $self->{instrument_model} = $instrument_model;
  }
  return $self->{instrument_model};
}

sub library_name{
  my ($self, $library_name) = @_;
  if($library_name){
    $self->{library_name} = $library_name;
  }
  return $self->{library_name};
}

sub run_name{
  my ($self, $run_name) = @_;
  if($run_name){
    $self->{run_name} = $run_name;
  }
  return $self->{run_name};
}

sub run_block_name{
  my ($self, $run_block_name) = @_;
  if($run_block_name){
    $self->{run_block_name} = $run_block_name;
  }
  return $self->{run_block_name};
}

sub paired_length{
  my ($self, $paired_length) = @_;
  if(defined($paired_length)){
    $self->{paired_length} = $paired_length;
  }
  return $self->{paired_length};
}

sub library_layout{
  my ($self, $library_layout) = @_;
  if($library_layout){
    $self->{library_layout} = $library_layout;
  }
  return $self->{library_layout};
}

sub status{
  my ($self, $status) = @_;
  if($status){
    $self->{status} = $status;
  }
  return $self->{status};
}

sub archive_base_count{
  my ($self, $archive_base_count) = @_;
  if(defined($archive_base_count)){
    $self->{archive_base_count} = $archive_base_count;
  }
  return $self->{archive_base_count};
}

sub archive_read_count{
  my ($self, $archive_read_count) = @_;
  if(defined($archive_read_count)){
    $self->{archive_read_count} = $archive_read_count;
  }
  return $self->{archive_read_count};
}

sub library_strategy{
  my ($self, $library_strategy) = @_;
  if(defined($library_strategy)){
    $self->{library_strategy} = $library_strategy;
  }
  return $self->{library_strategy};
}

sub fix_date{
  my ($self) = @_;
  if($self->submission_date &&
     $self->submission_date =~ /(\d\d)\-(\w\w\w)\-(\d\d\d\d)(.*)/){
    my $new = $3."-".$2."-".$1." ".$4;
    $self->submission_date($new);
  }
}

sub object_table_name{
  my ($self) = @_;
  return 'run_meta_info';
}

1;

