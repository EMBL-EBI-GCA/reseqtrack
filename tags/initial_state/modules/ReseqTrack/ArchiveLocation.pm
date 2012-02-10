
=pod

=head1 NAME

ReseqTrack::ArchiveLocation

=head1 SYNOPSIS

This is a container object for the archive_location table which provides human 
memorable names and machine useable path for the different archive_location ids

=cut

package ReseqTrack::ArchiveLocation;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

=head2 new

  Arg [1]   : ReseqTrack::Archive::Location
  Arg [2]   : string, path to root of location
  Arg [3]   : string, human name for location
  Function  : create ReseqTrack::ArchiveLocation object
  Returntype: ReseqTrack::ArchiveLocation
  Exceptions: throws if dbID isn't one of the accepted values
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my ( $location, $location_name ) =
      rearrange( [qw(LOCATION LOCATION_NAME)], @args );
    $self->location($location);
    $self->location_name($location_name);
    throw( "Can't have a dbID " . $self->dbID . " which isn't 1 or 2" )
      unless ( $self->dbID == 1 || $self->dbID == 2 );
    return $self;
}

=head2 Accessor methods

  Arg [1]   : ReseqTrack::ArchiveLocation
  Arg [2]   : string
  Function  : set internal hash value to passed in string
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub location {
    my ( $self, $location ) = @_;
    if ($location) {
        $self->{'location'} = $location;
    }
    return $self->{'location'};
}

sub location_name {
    my ( $self, $location_name ) = @_;
    if ($location_name) {
        $self->{'location_name'} = $location_name;
    }
    return $self->{'location_name'};
}


=head2 object_table_name

  Arg [1]   : ReseqTrack::ArchiveLocation
  Function  : return table name, archive_location
  Returntype: string
  Exceptions: n/a
  Example   : 

=cut



sub object_table_name {
    my ($self) = @_;
    return 'archive_location';
}

1;
