package ReseqTrack::Study;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($status, $md5, $type, $title ) =
      rearrange([qw(STATUS MD5 TYPE TITLE) ], @args);
  $self->status($status);
  $self->md5($md5);
  $self->type($type);
  $self->title($title);
  
  return $self;
}

sub status{
  my ($self, $arg) = @_; 
  
  if($arg){
    $self->{status} = $arg;
  }
  return $self->{status};
}

sub md5{
  my ($self, $arg) = @_;
  if($arg){
    $self->{md5} = $arg;
  }
  return $self->{md5};
}

sub type{
  my ($self, $arg) = @_;
  if($arg){
    $self->{type} = $arg;
  }
  return $self->{type};
}

sub title{
  my ($self, $arg) = @_;
  if($arg){
    $self->{title} = $arg;
  }
  return $self->{title};
}

1;
