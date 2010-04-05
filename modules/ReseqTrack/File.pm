
=pod

=head1 NAME

ReseqTrack::Event

=head1 SYNOPSIS

This is a container object for the file table. The file table describes the location, md5sum 
and size of the file alongside a type which defines a group of objects

=head1 Example

my $file = ReseqTrack::File->new(
      -name => $path,
      -type => $type,
      -size => $size,
      -host => $host,
        );

=cut

package ReseqTrack::File;

use strict;
use warnings;
use vars qw(@ISA);
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::HasHistory;

@ISA = qw(ReseqTrack::HasHistory);

=head2 new
  
  Arg [1]    : ReseqTrack::File
  Arg [2]    : string, name
  Arg [3]    : 32 bit string, md5
  Arg [4]    : string, type
  Arg [5]    : int, size of file
  Arg [6]    : string, directory path of file
  Arg [7]    : ReseqTrack::Host, object stating the host of the file
  Arg [8]    : int, binary value stating if the file is active in the project or not
  Arg [9]    : string, date time string when the file object was created
  Arg [10]   : string, date time string when the file object was updated
  Arg [11]   : int, database id for host object
  Function   : create a ReseqTrack::File object
  Returntype : ReseqTrack::File
  Exceptions : throws if no name, path, type or host are specified
  Example    : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my (
        $name,      $md5,     $type,    $size, $host,
        $withdrawn, $created, $updated, $host_id
      )
      = rearrange(
        [
            qw(NAME MD5 TYPE SIZE HOST WITHDRAWN CREATED UPDATED
              HOST_ID)
        ],
        @args
      );
    if ( $host_id && $self->adaptor ) {
        my $ha = $self->get_HostAdaptor;
        $host = $ha->fetch_by_dbID;
    }

    #ERROR CHECKING
    throw("Can't create ReseqTrack::File without a name") unless ($name);
    throw("Can't create ReseqTrack::File without a type") unless ($type);
    throw("Can't create ReseqTrack::File without a host object ")
      unless ( $host || $host->isa("ReseqTrack::Host") );
    ######
    $self->md5($md5);
    $self->type($type);
    $self->size($size);
    $self->name($name);
    $self->host($host);
    $self->withdrawn($withdrawn);
    $self->created($created);
    $self->updated($updated);
    #########

    if ( !$self->size && -e $self->full_path ) {
        my $size = -s $self->full_path;
        $self->size($size);
    }
    $self->withdrawn(0) unless ( defined( $self->withdrawn ) );
    return $self;
}

=head2 accessor methods

  Arg [1]   : ReseqTrack::File
  Arg [2]   : variable
  Function  : store variable in object
  Returntype: variable
  Exceptions: none
  Example   : 

=cut

sub type {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{type} = $arg;
    }
    return $self->{type};
}

sub size {
    my ( $self, $arg ) = @_;
    if ( defined($arg) ) {
        $self->{size} = $arg;
    }
    return $self->{size};
}

sub withdrawn {
    my ( $self, $arg ) = @_;
    if ( defined($arg) ) {
        $self->{withdrawn} = $arg;
    }
    return $self->{withdrawn};
}

sub created {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{created} = $arg;
    }
    return $self->{created};
}

sub updated {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{updated} = $arg;
    }
    return $self->{updated};
}

=head2 name

  Arg [1]   : ReseqTrack::File
  Arg [2]   : string, full file name and path
  Function  : stores path in object, will strip off terminating /
  Returntype: string, path
  Exceptions: none
  Example   : 

=cut

sub name {
    my ( $self, $arg ) = @_;
    if ($arg) {

        #trying to prevent database / and // issues
        $arg =~ s/\/$//;
        $arg =~ s/\/\//\//;
        $self->{name} = $arg;
    }
    return $self->{name};
}

=head2 md5

  Arg [1]   : ReseqTrack::File
  Arg [2]   : string, md5sum value should be 32 characters long
  Function  : storing md5 value in object
  Returntype: string
  Exceptions: throws if string isn't 32 characters long
  Example   : 

=cut

sub md5 {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{md5} = $arg;
    }
    if ( $self->{md5} && length( $self->{md5} ) != 32 ) {
        throw(  "Have give ReseqTrack::File::md5 a string which is "
              . length( $self->{md5} )
              . " not 32 characters long" );
    }
    return $self->{md5};
}

=head2 host

  Arg [1]   : ReseqTrack::File
  Arg [2]   : ReseqTrack::Host
  Function  : storing host object in current object
  Returntype: ReseqTrack::Host
  Exceptions: throw if object passed in isn't a ReseqTrack::Host
  Example   : 

=cut

sub host {
    my ( $self, $host ) = @_;
    if ($host) {
        throw("ReseqTrack::File::host must be a ReseqTrack::Host object")
          unless ( $host->isa("ReseqTrack::Host") );
        $self->{host} = $host;
    }
    return $self->{host};
}

=head2 filename/path

  Arg [1]   : ReseqTrack::File
  Function  : retrun the filename or directory of the file path
  Returntype: string
  Exceptions: none
  Example   : my $filename = $file->filename;

=cut

sub filename {
    my ($self) = @_;
    if ( !$self->{filename} ) {
        my $filename = basename( $self->name );
        $self->{filename} = $filename;
    }
    return $self->{filename};
}

sub path {
    my ($self) = @_;
    if ( !$self->{dir} ) {
        my $dir = dirname( $self->name );
        $self->{dir} = $dir;
    }
    return $self->{dir};
}

=head2 full_path

  Arg [1]   : ReseqTrack::File
  Function  : return $self->path
  Returntype: string
  Exceptions: n/a
  Example   : 

=cut

sub full_path {
    my ($self) = @_;
    return $self->name;
}

=head2 object_table_name

  Arg [1]   : ReseqTrack::Event
  Function  : return table name for object, event
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub object_table_name {
    my ($self) = @_;
    return 'file';
}

1;
