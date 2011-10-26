=pod

=head1 NAME

ReseqTrack::Archive

=head1 SYNOPSIS

This is a container object for the archive table. The archive table is a system
to interact with the fire archive system and allow data to be moved to the replicated
archive system that the systems team runs

=head1 Example

my $archive = ReseqTrack::Archive->new
      (
       -file => $file,
       -archive_action => $archive_action,
       -archive_location => $archive_location,
      );

=cut

package ReseqTrack::Archive;

use strict;
use warnings;
use vars qw(@ISA);

use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

=head2 new

  Arg [1]   : ReseqTrack::Archive
  Arg [2]   : string, file name
  Arg [3]   : int, file dbID
  Arg [4]   : string , md5 checksum
  Arg [5]   : string, relative path to file from archive_location
  Arg [6]   : string, the name of the volume is the archive location has one
  Arg [7]   : string, timestamp for when created
  Arg [8]   : string, timestamp for when updated
  Arg [9]   : int, archive action dbID
  Arg [10]  : ReseqTrack::File
  Arg [11]  : ReseqTrack::ArchiveAction
  Arg [12]  : ReseqTrack::ArchiveLocation
  Arg [13]  : string, new name
  Arg [14]  : string, new relative path
  Arg [15]  : int, priority, must be between 1-100
  Arg [16]  : int, fire action id (used by the fire system)
  Arg [17]  : int, fire exit code, if there is a problem
  Arg [18]  : string, fire exit reason, explaination of problem
  Function  : create a ReseqTrack::Archive object
  Returntype: ReseqTrack::Archive
  Exceptions: throws if no archive location is set
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my (
        $name,           $file_id,           $md5,
        $relative_path,  $volume_name,       $created,
        $updated,        $archive_action_id, $archive_location_id,
        $file,           $archive_action,    $archive_location,
        $new_name,       $new_relative_path, $priority,
        $fire_action_id, $fire_exit_code,    $fire_exit_reason,
        $size
      )
      = rearrange(
        [
            'NAME',                'FILE_ID',
            'MD5',                 'RELATIVE_PATH',
            'VOLUME_NAME',         'CREATED',
            'UPDATED',             'ARCHIVE_ACTION_ID',
            'ARCHIVE_LOCATION_ID', 'FILE',
            'ARCHIVE_ACTION',      'ARCHIVE_LOCATION',
            'NEW_NAME',            'NEW_RELATIVE_PATH',
            'PRIORITY',            'FIRE_ACTION_ID',
            'FILE_EXIT_CODE',      'FILE_EXIT_REASON',
            'SIZE'
        ],
        @args
      );

    $self->name($name);
    $self->file_id($file_id);
    $self->md5($md5);
    $self->relative_path($relative_path);
    $self->volume_name($volume_name);
    $self->created($created);
    $self->updated($updated);
    $self->archive_action_id($archive_action_id);
    $self->archive_location_id($archive_location_id);
    $self->archive_action($archive_action);
    $self->archive_location($archive_location);
    $self->new_name($new_name);
    $self->new_name($name) if ( !$self->new_name );
    $self->new_relative_path($new_relative_path);
    $self->priority($priority);
    $self->fire_action_id($fire_action_id);
    $self->fire_exit_code($fire_exit_code);
    $self->fire_exit_reason($fire_exit_reason);
    throw("Can't archive an object without at least an archive location id")
      unless ( $self->archive_location_id );
    $self->file($file);
    $self->size($size);
    $self->priority(50) unless ( $self->priority );
    return $self;
}

=head2 file

  Arg [1]   : ReseqTrack::Archive
  Arg [2]   : ReseqTrack::File
  Function  : Uses a ReseqTrack File object to define other variables required
  by the archive system or can use the specified file id to retrieve the
  appropriate file object
  Returntype: ReseqTrack::File
  Exceptions: 
  Example   : 

=cut

sub file {
    my ( $self, $object ) = @_;
    if ($object) {
        throw( $object . " must be a ReseqTrack::File object" )
          unless ( $object->isa('ReseqTrack::File') );
        $self->{'file'} = $object;
    }
    if ( $self->{'file'}
        && ( !$self->name && !$self->md5 && !$self->relative_path ) )
    {
        $self->name( $self->{'file'}->filename );
        $self->md5( $self->{'file'}->md5 );
        $self->file_id( $self->{'file'}->dbID );
        $self->size( $self->{'file'}->size );
        my $relative_path = $self->{'file'}->full_path;
        $relative_path = dirname($relative_path);
        my $root = $self->archive_location->location;
        $relative_path =~ s/$root//;
        my $vol_name;

        if ( $relative_path =~ /vol\d+/ ) {
            $relative_path =~ /(vol\d+)\/(.+)/i;
            $vol_name      = $1;
            $relative_path = $2;
            if ( !$vol_name || !$relative_path ) {
                throw( "Couldn't parse relative_path from "
                      . $self->{'file'}->full_path );
            }
        }
        $self->volume_name($vol_name);
        $self->relative_path($relative_path);
    }
    if ( !$self->{'file'} && $self->adaptor && $self->file_id ) {
        $self->{'file'} =
          $self->adaptor->db->get_FileAdaptor->fetch_by_dbID( $self->file_id );
    }
    return $self->{'file'};
}

=head2 archive_action/location

  Arg [1]   : ReseqTrack::Archive
  Arg [2]   : ReseqTrack::ArchiveAction/ArchiveLocation
  Function  : acessor method for ArchiveAction or ArchiveLocation objects
  will use defined id to fetch given objects if not already present
  Returntype: ReseqTrack::ArchiveAction/ArchiveLocation
  Exceptions: throws if not given the correct object type
  Example   : 

=cut

sub archive_action {
    my ( $self, $object ) = @_;
    if ($object) {
        throw( $object . " must be a ReseqTrack::ArchiveAction object" )
          unless ( $object->isa('ReseqTrack::ArchiveAction') );
        $self->{'archive_action'} = $object;
    }
    if ( !$self->archive_action_id && $self->{archive_action} ) {
        $self->archive_action_id( $self->{archive_action}->dbID );
    }
    if (  !$self->{archive_action}
        && $self->archive_action_id
        && $self->adaptor )
    {
        $self->{archive_action} =
          $self->adaptor->db->get_ArchiveActionAdaptor->fetch_by_dbID(
            $self->archive_action_id );
    }
    return $self->{'archive_action'};
}

sub archive_location {
    my ( $self, $object ) = @_;
    if ($object) {
        throw( $object . " must be a ReseqTrack::ArchiveLocation object" )
          unless ( $object->isa('ReseqTrack::ArchiveLocation') );
        $self->{archive_location} = $object;
    }
    if ( !$self->archive_location_id && $self->{archive_location} ) {
        $self->archive_location_id( $self->{archive_location}->dbID );
    }
    if (  !$self->{archive_location}
        && $self->archive_location_id
        && $self->adaptor )
    {
        $self->{archive_location} =
          $self->adaptor->db->get_ArchiveLocationAdaptor->fetch_by_dbID(
            $self->archive_location_id );
    }
    return $self->{'archive_location'};
}

=head2 Accessor methods

  Arg [1]   : ReseqTrack::Archive
  Arg [2]   : string/int
  Function  : set/return given variable
  Returntype: string/int
  Exceptions: 
  Example   : 

=cut

sub created {
    my ( $self, $time ) = @_;
    if ($time) {
        $self->{'created'} = $time;
    }
    return $self->{'created'};
}

sub updated {
    my ( $self, $time ) = @_;
    if ($time) {
        $self->{'updated'} = $time;
    }
    return $self->{'updated'};
}

sub name {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'name'} = $value;
    }
    return $self->{'name'};
}

sub new_name {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'new_name'} = $value;
    }
    return $self->{'new_name'};
}

sub file_id {

    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'file_id'} = $value;
    }
    return $self->{'file_id'};
}

sub md5 {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'md5'} = $value;
    }
    return $self->{'md5'};
}

sub size {
    my ( $self, $value ) = @_;
    if ( defined ($value) ) {
        $self->{'size'} = $value;
    }
    if ( !$self->{'size'} ) {
        my $full_path = $self->full_path;
        $self->{'size'} = -s $full_path;
    }
    return $self->{'size'};
}

sub priority {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'priority'} = $value;
    }
    return $self->{'priority'};
}

sub relative_path {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'relative_path'} = $value;
    }
    return $self->{'relative_path'};
}

sub new_relative_path {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'new_relative_path'} = $value;
    }
    return $self->{'new_relative_path'};
}



sub volume_name {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'volume_name'} = $value;
    }
    return $self->{'volume_name'};
}

sub archive_action_id {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'archive_action_id'} = $value;
    }
    return $self->{'archive_action_id'};
}

sub archive_location_id {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'archive_location_id'} = $value;
    }
    return $self->{'archive_location_id'};
}

sub fire_action_id {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'fire_action_id'} = $value;
    }
    return $self->{'fire_action_id'};
}

sub fire_exit_code {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'fire_exit_code'} = $value;
    }
    return $self->{'fire_exit_code'};
}

sub fire_exit_reason {
    my ( $self, $value ) = @_;
    if ($value) {
        $self->{'fire_exit_reason'} = $value;
    }
    return $self->{'fire_exit_reason'};
}

=head2 full_path

  Arg [1]   : ReseqTrack::Archive
  Function  : create the full path to the given file based on the information
  held by the object
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub full_path {
    my ($self) = @_;
    my $path = $self->archive_location->location;
    $path .= "/" unless ( $path =~ /\/$/ );
    $path .= $self->volume_name if ( $self->volume_name );
    $path .= "/" unless ( $path =~ /\/$/ );
    $path .= $self->relative_path;
    $path .= "/" unless ( $path =~ /\/$/ );
    $path .= $self->name;
    $path =~ s/\/\//\//g;
    return $path;
}

=head2 object_table_name

  Arg [1]   : ReseqTrack::Archive
  Function  : return table name for object, archive
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub table_name {
    my ($self) = @_;
    return "archive";
}

1;

