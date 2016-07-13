package ReseqTrack::Sample;

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

  my (
    $status,          $md5,             $center_name,
    $sample_alias,    $tax_id,          $scientific_name,
    $common_name,     $anonymized_name, $individual_name,
    $sample_title,    $source_id,       $submission_id,
    $submission_date, $biosample_id,    $biosample_authority
    )
    = rearrange(
    [
      qw(STATUS MD5 CENTER_NAME
        SAMPLE_ALIAS TAX_ID SCIENTIFIC_NAME COMMON_NAME
        ANONYMIZED_NAME INDIVIDUAL_NAME SAMPLE_TITLE SOURCE_ID SUBMISSION_ID SUBMISSION_DATE BIOSAMPLE_ID
        BIOSAMPLE_AUTHORITY)
    ],
    @args
    );

  $self->status($status);
  $self->md5($md5);
  $self->center_name($center_name);
  $self->sample_alias($sample_alias);
  $self->tax_id($tax_id);
  $self->scientific_name($scientific_name);
  $self->common_name($common_name);
  $self->anonymized_name($anonymized_name);
  $self->individual_name($individual_name);
  $self->sample_title($sample_title);
  $self->source_id($source_id);
  $self->submission_id($submission_id);
  $self->submission_date($submission_date);
  $self->biosample_id($biosample_id);
  $self->biosample_authority($biosample_authority);

  return $self;
}

sub name {
  my ($self) = @_;
  return $self->source_id();
}

sub source_id {
  my ( $self, $arg ) = @_;
  return $self->sample_source_id($arg);
}

sub sample_source_id {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{sample_source_id} = $arg;
  }
  return $self->{sample_source_id};
}

sub status {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{status} = $arg;
  }
  return $self->{status};
}

sub biosample_id {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{biosample_id} = $arg;
  }
  return $self->{biosample_id};
}

sub biosample_authority {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{biosample_authority} = $arg;
  }
  return $self->{biosample_authority};
}

sub md5 {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{md5} = $arg;
  }
  return $self->{md5};
}

sub center_name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{center_name} = $arg;
  }
  return $self->{center_name};
}

sub sample_alias {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{sample_alias} = $arg;
  }
  return $self->{sample_alias};
}

sub tax_id {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{tax_id} = $arg;
  }
  return $self->{tax_id};
}

sub scientific_name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{scientific_name} = $arg;
  }
  return $self->{scientific_name};
}

sub common_name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{common_name} = $arg;
  }
  return $self->{common_name};
}

sub anonymized_name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{anonymized_name} = $arg;
  }
  return $self->{anonymized_name};
}

sub individual_name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{individual_name} = $arg;
  }
  return $self->{individual_name};
}

sub sample_title {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{sample_title} = $arg;
  }
  return $self->{sample_title};
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
  return "sample";
}

1;
