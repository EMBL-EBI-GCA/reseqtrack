package ReseqTrack::Study;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::HasHistory;

@ISA = qw(ReseqTrack::HasHistory);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $status, $md5, $type, $title, $alias, $source_id, $submission_id,
    $submission_date )
    = rearrange(
    [qw(STATUS MD5 TYPE TITLE STUDY_ALIAS SOURCE_ID SUBMISSION_ID SUBMISSION_DATE)],
    @args );
  $self->status($status);
  $self->md5($md5);
  $self->type($type);
  $self->title($title);
  $self->study_alias($alias);
  $self->source_id($source_id);
  $self->submission_id($submission_id);
  $self->submission_date($submission_date);

  return $self;
}

sub name {
  my ($self) = @_;
  return $self->source_id();
}

sub source_id {
  my ( $self, $arg ) = @_;
  return $self->study_source_id($arg);
}

sub study_source_id {
	my ($self,$arg) =@_;
	
	if ($arg) {
    $self->{study_source_id} = $arg;
  }
  return $self->{study_source_id};
}


sub status {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{status} = $arg;
  }
  return $self->{status};
}

sub md5 {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{md5} = $arg;
  }
  return $self->{md5};
}

sub type {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{type} = $arg;
  }
  return $self->{type};
}

sub title {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{title} = $arg;
  }
  return $self->{title};
}
sub study_alias {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{study_alias} = $arg;
  }
  return $self->{study_alias};
}

sub submission_id {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{submission_id} = $arg;
  }
  return $self->{submission_id};
}

sub submission_date {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{submission_date} = $arg;
  }
  return $self->{submission_date};
}

sub object_table_name {
  return "study";
}

1;
