package ReseqTrack::InputString;

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
  
  my ($name, $type ) = rearrange([qw( NAME TYPE ) ], @args);

  $self->name($name);
  $self->type($type);
  throw("Must have a name for a ReseqTrack::InputString") unless($name);
  return $self;
}

#name
#type

sub name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{name} = $arg;
  }
  return $self->{name};
}

sub type{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{type} = $arg;
  }
  return $self->{type};
}


1;
