package ReseqTrack::Run;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Base;

@ISA = qw(ReseqTrack::HasHistory);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my (
    $experiment_id,       $run_alias,        $status,
    $md5,                 $center_name,      $run_center_name,
    $instrument_platform, $instrument_model, $source_id,
    $sample_id,           $submission_id, $submission_date
    )
    = rearrange(
    [
      qw(
        EXPERIMENT_ID
        RUN_ALIAS
        STATUS
        MD5
        CENTER_NAME
        RUN_CENTER_NAME
        INSTRUMENT_PLATFORM
        INSTRUMENT_MODEL
        SOURCE_ID
        SAMPLE_ID
        SUBMISSION_ID
        SUBMISSION_DATE)
    ],
    @args
    );

  $self->experiment_id($experiment_id);
  $self->run_alias($run_alias);
  $self->status($status);
  $self->md5($md5);
  $self->center_name($center_name);
  $self->run_center_name($run_center_name);
  $self->instrument_platform($instrument_platform);
  $self->instrument_model($instrument_model);
  $self->source_id($source_id);
  $self->sample_id($sample_id);
  $self->submission_id($submission_id);
  $self->submission_date($submission_date);
  return $self;
}


sub experiment {
  my ($self) = @_;
  
  if ( $self->adaptor && $self->experiment_id ) {
      my $ea = $self->adaptor->db->get_ExperimentAdaptor;
      return $ea->fetch_by_dbID($self->experiment_id);
  }
}

sub sample {
  my ($self) = @_;
  
  if ( $self->adaptor && $self->sample_id ) {
      my $sa = $self->adaptor->db->get_SampleAdaptor;
      return $sa->fetch_by_dbID($self->sample_id);
  }
}

sub name {
  my ($self) = @_;
  return $self->source_id();
}

sub source_id {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{source_id} = $arg;
  }
  return $self->{source_id};
}

sub submission_id{
  my ($self, $arg) = @_;
  if($arg){
    $self->{submission_id} = $arg;
  }
  return $self->{submission_id};
}

sub submission_date{
  my ($self, $arg) = @_;
  if($arg){
    $self->{submission_date} = $arg;
  }
  return $self->{submission_date};
}

sub experiment_id {
  my ( $self, $arg ) = @_;
  if ($arg) { $self->{experiment_id} = $arg; }
  return $self->{experiment_id};
}

sub run_alias {
  my ( $self, $arg ) = @_;
  if ($arg) { $self->{run_alias} = $arg; }
  return $self->{run_alias};
}

sub status {
  my ( $self, $arg ) = @_;
  if ($arg) { $self->{status} = $arg; }
  return $self->{status};
}

sub md5 {
  my ( $self, $arg ) = @_;
  if ($arg) { $self->{md5} = $arg; }
  return $self->{md5};
}

sub center_name {
  my ( $self, $arg ) = @_;
  if ($arg) { $self->{center_name} = $arg; }
  return $self->{center_name};
}

sub run_center_name {
  my ( $self, $arg ) = @_;
  if ($arg) { $self->{run_center_name} = $arg; }
  return $self->{run_center_name};
}

sub instrument_platform {
  my ( $self, $arg ) = @_;
  if ($arg) { $self->{instrument_platform} = $arg; }
  return $self->{instrument_platform};
}

sub instrument_model {
  my ( $self, $arg ) = @_;
  if ($arg) { $self->{instrument_model} = $arg; }
  return $self->{instrument_model};
}

sub sample_id {
  my ( $self, $arg ) = @_;
  if ($arg) { $self->{sample_id} = $arg; }
  return $self->{sample_id};
}

sub object_table_name {
  return "run";
}
