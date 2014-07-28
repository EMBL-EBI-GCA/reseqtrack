package ReseqTrack::Tools::Counts;

# Simple value class to hold the base and read counts for a Fastq file, before and after QA
# Holly Zheng Bradley, Oct.  2009

use strict;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);

sub new{
  my ($class, @args) = @_;
  my $self = bless( {}, $class);

  my ($file, $base_cnt, $read_cnt, $filt_base_cnt, $filt_read_cnt) =
      rearrange([qw(FILE BASE_CNT READ_CNT FILT_BASE_CNT FILT_READ_CNT)], @args);

  ## initiating the object
  $self->file($file);
  $self->base_cnt($base_cnt);
  $self->read_cnt($read_cnt);
  $self->filt_base_cnt($filt_base_cnt);
  $self->filt_read_cnt($filt_read_cnt);
  return $self;
}

sub file{
  my ($self, $arg) = @_;
  $self->{file} = $arg if $arg;
  return $self->{file};
}

sub base_cnt{
  my ($self, $arg) = @_;
  $self->{base_cnt} = $arg if $arg;
  return $self->{base_cnt};
}

sub read_cnt{
  my ($self, $arg) = @_;
  $self->{read_cnt} = $arg if $arg;
  return $self->{read_cnt};
}

sub filt_base_cnt{
  my ($self, $arg) = @_;
  $self->{filt_base_cnt} = $arg if $arg;
  return $self->{filt_base_cnt};
}

sub filt_read_cnt{
  my ($self, $arg) = @_;
  $self->{filt_read_cnt} = $arg if $arg;
  return $self->{filt_read_cnt};
}

1;
