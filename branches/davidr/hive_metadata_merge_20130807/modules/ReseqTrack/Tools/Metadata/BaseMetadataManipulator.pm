package ReseqTrack::Tools::Metadata::BaseMetadataManipulator;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);

sub new {
  my ($class, @args) = @_;
  my $self ={};
  bless $self,$class;
  my ($era_db, $reseq_db,) = rearrange([qw(ERA_DB RESEQ_DB)], @args);

  $self->era_db($era_db);
  $self->reseq_db($reseq_db);

  return $self;
}

sub era_db {
  my ($self, $arg) = @_; 
  
  if($arg){
    $self->{era_db} = $arg;
  }
  return $self->{era_db};
}
sub reseq_db {
  my ($self, $arg) = @_; 
  
  if($arg){
    $self->{reseq_db} = $arg;
  }
  return $self->{reseq_db};
}
sub manipulate {
  my ($self,$object,$last_copy) = @_;
  
  if ($object->isa('ReseqTrack::Study')) {
    return $self->manipulate_study($object,$last_copy);
  }
  if ($object->isa('ReseqTrack::Sample')) {
    return $self->manipulate_sample($object,$last_copy);
  }
  if ($object->isa('ReseqTrack::Experiment')) {
    return $self->manipulate_experiment($object,$last_copy);
  }
  if ($object->isa('ReseqTrack::Run')) {
    return $self->manipulate_run($object,$last_copy);
  }
  
}

sub manipulate_study{}
sub manipulate_sample{}
sub manipulate_experiment{}
sub manipulate_run{}


1;