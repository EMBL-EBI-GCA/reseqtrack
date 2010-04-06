
=pod

=head1 NAME

ReseqTrack::Host

=head1 SYNOPSIS

A container to hold the data for the host table in the ReseqTrack Database

The host table describes the location of the host machine the file lives on and
if that machine is local or remote

=head1 EXAMPLE

my $host = ReseqTrack::Host->new(
               -name => '1000genomes.ebi.ac.uk', 
               -remote => 0,
             );

=cut

package ReseqTrack::Host;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);


=head2 new

  Arg [1]   : ReseqTrack::Host
  Arg [2]   : string, name of host
  Arg [3]   : binary 0/1 set to 0 if the host is local to where the database is
  Function  : create ReseqTrack object
  Returntype: ReseqTrack::Host
  Exceptions: throws if no name is defined
  Example   : 

=cut



sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $name, $remote ) = rearrange( [qw( NAME REMOTE )], @args );

    #setting default
    $self->remote(0);
    ##########
    $self->name($name);
    $self->remote($remote);
    throw("Must have a name for a ReseqTrack::Host") unless ($name);
    return $self;
}



=head2 Accessor methods

  Arg [1]   : Reseq::Host
  Arg [2]   : mostly strings
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

sub remote {
    my ( $self, $arg ) = @_;
    if ( defined($arg) ) {
        $self->{remote} = $arg;
    }
    return $self->{remote};
}

1;
