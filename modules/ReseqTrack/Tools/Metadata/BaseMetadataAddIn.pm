package ReseqTrack::Tools::Metadata::BaseMetadataAddIn;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);

sub new {
  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;
  my ( $era_db, $reseq_db, $log_fh ) =
    rearrange( [qw(ERA_DB RESEQ_DB LOG_FH)], @args );

  $self->era_db($era_db);
  $self->reseq_db($reseq_db);
  $self->log_fh($log_fh);

  return $self;
}

sub era_db {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{era_db} = $arg;
  }
  return $self->{era_db};
}

sub reseq_db {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{reseq_db} = $arg;
  }
  return $self->{reseq_db};
}

sub log_fh {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{log_fh} = $arg;
  }
  return $self->{log_fh};
}

sub check {
  my ( $self, $object, $last_copy ) = @_;

  if ( $object->isa('ReseqTrack::Study') ) {
    return $self->check_study( $object, $last_copy );
  }
  if ( $object->isa('ReseqTrack::Sample') ) {
    return $self->check_sample( $object, $last_copy );
  }
  if ( $object->isa('ReseqTrack::Experiment') ) {
    return $self->check_experiment( $object, $last_copy );
  }
  if ( $object->isa('ReseqTrack::Run') ) {
    return $self->check_run( $object, $last_copy );
  }

}

sub check_study      { }
sub check_sample     { }
sub check_experiment { }
sub check_run        { }
sub report           { }

1;
