package ReseqTrack::Host;

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
  
  my ($name, $remote ) = rearrange([qw( NAME REMOTE ) ], @args);

  #setting default
  $self->remote(0);
  ##########
  $self->name($name);
  $self->remote($remote);
  throw("Must have a name for a ReseqTrack::Host") unless($name);
  return $self;
}

#name
#remote

sub name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{name} = $arg;
  }
  return $self->{name};
}

sub remote{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{remote} = $arg;
  }
  return $self->{remote};
}


1;
