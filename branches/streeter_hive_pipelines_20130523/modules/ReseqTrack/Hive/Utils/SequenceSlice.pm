

package ReseqTrack::Hive::Utils::SequenceSlice;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use POSIX qw(ceil);

sub new {
    my ( $class, @args ) = @_;
    my $self = {};
    bless $self, $class;

    my ( $SQ_name, $SQ_length, $start, $end)
      = rearrange( [ qw(SQ_NAME SQ_LENGTH START END ) ],
        @args
      );

    throw("SequenceSlice object must know the SQ_name") if (!$SQ_name);
    throw("SequenceSlice object must know the SQ_length") if (!$SQ_length);
    
    $self->SQ_name($SQ_name);
    $self->SQ_length($SQ_length);
    $self->start($start || 1);
    $self->end($end || $SQ_length);

    return $self;
}

sub is_whole_SQ {
    my ( $self, $arg ) = @_;
    if ($arg) {
      $self->start(1);
      $self->end($self->SQ_length);
    }
    return ($self->start == 1 && $self->end == $self->SQ_length) ? 1 : 0;
}

sub SQ_name {
    my ( $self, $arg ) = @_;
    if ($arg) {
      $self->{SQ_name} = $arg;
    }
    return $self->{SQ_name};
}

sub SQ_length {
    my ( $self, $arg ) = @_;
    if ($arg) {
      $self->{SQ_length} = $arg;
    }
    return $self->{SQ_length};
}

sub start {
    my ( $self, $arg ) = @_;
    if (defined $arg) {
        $self->{start} = $arg;
    }
    return $self->{start};
}

sub end {
    my ( $self, $arg ) = @_;
    if (defined $arg) {
        $self->{end} = $arg;
    }
    return $self->{end};
}

sub length {
    my ( $self ) = @_;
    return $self->end - $self->start +1;
}

sub extend_5prime {
  my ($self, $extra_length) = @_;

  my $start = $self->start - $extra_length;
  $start = $start < 1 ? 1 : $start;
  $self->start($start);
}

sub extend_3prime {
  my ($self, $extra_length) = @_;

  my $end = $self->end + $extra_length;
  $end = $end > $self->SQ_length ? $self->SQ_length : $end;
  $self->end($end);
}

sub extend {
  my ($self, $extra_length) = @_;
  $self->extend_3prime($extra_length);
  $self->extend_5prime($extra_length);
}

sub split {
  my ($self, $max_length) = @_;
  my @slices;
  my $num_slices = ceil( ($self->length ) / $max_length );
  return [$self] if $num_slices == 1;
  my $slice_length = ceil( ($self->length ) / $num_slices);
  foreach my $i (0..$num_slices-1) {
    my $slice_start = $self->start + $i*$slice_length;
    my $slice_end = ($i == $num_slices-1) ? $self->end : $slice_start + $slice_length -1;
    my $slice = ReseqTrack::Hive::Utils::SequenceSlice->new(
        -SQ_name => $self->SQ_name,
        -SQ_length => $self->SQ_length,
        -start => $slice_start,
        -end => $slice_end,
        );
    push(@slices, $slice);
  }
  return \@slices;
}

1;
