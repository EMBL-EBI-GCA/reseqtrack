
=pod

=head1 NAME

ReseqTrack::InputString

=head1 SYNOPSIS

A container to hold the data for the input_string table in the ReseqTrack Database

The input string tables allows strings which can not be defined in either the 
collection, run meta info or file table to be used as input for event objects. The
name specifies the string and the type which set of inputs it belongs with

=head1 EXAMPLE

my $host = ReseqTrack::InputString->new(
               -name => 'chromosome1', 
               -type => 'CHROMOSOMES',
             );

=cut

package ReseqTrack::InputString;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

=head2 new

  Arg [1]   : ReseqTrack::Input string
  Arg [2]   : string, input string
  Arg [3]   : string, type string
  Function  : create ReseqTrack InputString object
  Returntype: ReseqTrack::InputString
  Exceptions: throws if no name is defined
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $name, $type ) = rearrange( [qw( NAME TYPE )], @args );

    $self->name($name);
    $self->type($type);
    throw("Must have a name for a ReseqTrack::InputString") unless ($name);
    return $self;
}

=head2 Accessor methods

  Arg [1]   : Reseq::InputString
  Arg [2]   : string
  Function  : set variable in object
  Returntype: return variable
  Exceptions: n/a
  Example   : 

=cut

sub name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{name} = $arg;
  }
  return $self->{name};
}

sub type {
  my ( $self, $arg ) = @_;
  if ( defined($arg) ) {
    $self->{type} = $arg;
  }
  return $self->{type};
}

1;
