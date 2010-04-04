
=pod

=head1 NAME

ReseqTrack::EventComplete

=head1 SYNOPSIS

This object is used to represent finished events and link those to specific inputs

=head1 Example

 my $completed_string = ReseqTrack::EventComplete->new
       (
        -event => $event,
        -other_name => $job->input_string,
        -sucess => 0,
        -adaptor => $ca,
       );

=cut

package ReseqTrack::EventComplete;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Base;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::Base);

=head2 new

  Arg [1]   : ReseqTrack::EventComplete
  Arg [2]   : ReseqTrack::Event
  Arg [3]   : int, event dbID
  Arg [4]   : ReseqTrack::Base, this can be another object has a dbID associated with it
  Arg [5]   : int, other dbID
  Arg [6]   : string, type this should be the type associated with the event
  Arg [7]   : string, table name this should be the table name associated with the input
  Arg [8]   : binary, 0/1 to indicate sucess (not yet properly implemented)
  Arg [9]   : time, the time the object was created
  Arg [10]  : string, other name this should be the name associated with the other object
  Function  : create EventComplete object
  Returntype: ReseqTrack::EventComplete
  Exceptions: throws if table name doesn't match expectation and if dbIDs have
  been passed without an adaptor
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my (
        $event,      $event_id, $other, $other_id, $type,
        $table_name, $success,  $time,  $other_name
      )
      = rearrange(
        [
            qw(EVENT EVENT_ID OTHER OTHER_ID TYPE TABLE_NAME
              SUCCESS TIME OTHER_NAME)
        ],
        @args
      );

    #error checking
    $success = 0 unless ( defined($success) );
    $table_name = $event->table_name if ( $event && !$table_name );
    $type = $event->type if ( !$type );
    throw(  "Must pass EventComplete a table name and it must be file "
          . "or collection or run_meta_info" )
      if (
        !$table_name
        || (   $table_name ne 'file'
            && $table_name ne 'collection'
            && $table_name ne 'run_meta_info'
            && $table_name ne 'input_string' )
      );
    if ( ( $event_id || $other_id ) && !( $self->adaptor ) ) {
        throw(  "EventComplete must have an adaptor if give event_id or "
              . $table_name
              . "_id" );
    }
    throw("EventComplete must be given a type") if ( !$type );
    throw("Can't create an EventComplete object without an adaptor")
      if ( !$self->adaptor );
    ######
    $self->table_name($table_name);
    $self->event_id($event_id);
    $self->event($event);
    $self->other($other);
    $self->other_id($other_id);
    $self->other_name($other_name);
    $self->type($type);
    $self->success($success);
    $self->time($time);
    return $self;
}

=head2 Accessor methods

  Arg [1]   : ReseqTrack::Archive
  Arg [2]   : string/int
  Function  : set/return given variable
  Returntype: string/int
  Exceptions: 
  Example   : 

=cut

sub table_name {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'table_name'} = $arg;
    }
    return $self->{'table_name'};
}

sub type {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'type'} = $arg;
    }
    return $self->{'type'};
}

sub success {
    my ( $self, $arg ) = @_;
    if ( defined($arg) ) {
        $self->{'success'} = $arg;
    }
    return $self->{'success'};
}

sub time {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'time'} = $arg;
    }
    return $self->{'time'};
}

=head2 event

  Arg [1]   : ReseqTrack::EventComplete
  Arg [2]   : ReseqTrack::Event
  Function  : store and return ReseqTrack event object
  Returntype: ReseqTrack::Event
  Exceptions: throw is not passed a ReseqTrack::Event
  Example   : 

=cut

sub event {
    my ( $self, $arg ) = @_;
    if ($arg) {
        throw( "Must pass ReseqTrack::EventComplete::event an event object not "
              . $arg )
          unless ( $arg->isa("ReseqTrack::Event") );
        $self->{'event'} = $arg;
    }
    return $self->{'event'};
}

=head2 event_id

  Arg [1]   : ReseqTrack::EventComplete
  Arg [2]   : int, dbID for Event object
  Function  : store event object dbID, if an adaptor is defined but an event object isn't
  Returntype: int, dbID for Event object
  Exceptions: throw if it no event object is defined and the dbID can't be used to find one
  Example   : 

=cut

sub event_id {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'event_id'} = $arg;
    }
    if ( $self->adaptor ) {
        if ( !$self->event ) {
            my $ea    = $self->adaptor->db->get_EventAdaptor;
            my $event = $ea->fetch_by_dbID( $self->{'event_id'} );
            throw( "Failed to fetch an event with dbID " . $self->{'event_id'} )
              if ( $self->{'event_id'} && !$event );
            $self->event($event);
        }
    }
    if ( $self->event && !$self->{'event_id'} ) {
        $self->{'event_id'} = $self->event->dbID;
    }
    return $self->{'event_id'};
}

=head2 other

  Arg [1]   : ReseqTrack::EventComplete
  Arg [2]   : ReseqTrack::Base
  Function  : holder for the object which provides the input string
  Returntype: ReseqTrack::Base
  Exceptions: 
  Example   : 

=cut

sub other {
    my ( $self, $other ) = @_;
    if ($other) {
        $self->{other}      = $other;
        $self->{other_id}   = $other->dbID;
        $self->{other_name} = $other->name;
    }
    return $self->{other};
}

=head2 other_id

  Arg [1]   : ReseqTrack::Collection
  Arg [2]   : int, dbID of other object
  Function  : hold other dbID and use dbID to fetch an other object if other object isn't 
  already defined
  Returntype: int
  Exceptions: throw if it can't fetch the object using the given dbID
  Example   : 

=cut

sub other_id {
    my ( $self, $id ) = @_;
    if ($id) {
        $self->{other_id} = $id;
        if ( $self->adaptor && !$self->other ) {
            my $oa    = $self->get_other_adaptor;
            my $other = $oa->fetch_by_dbID($id);
            warning(
                "Failed to fetch " . $self->table_name . " using id " . $id )
              if ( !$other );
            $self->{other} = $other;
            $self->{other_name} = $other->name if ($other);
        }
    }
    return $self->{other_id};
}

=head2 other_name

  Arg [1]   : ReseqTrack::EventComplete
  Arg [2]   : string, name of other object
  Function  : stores name or uses object to defined object name
  Returntype: string
  Exceptions: throws if no other object is defined so the name can't be fetched
  Example   : 

=cut

sub other_name {
    my ( $self, $name ) = @_;
    if ($name) {
        $self->{other_name} = $name;
        if ( $self->adaptor && !$self->other ) {
            if ( !$self->event ) {
                throw("Don't have an event object, not sure what to do");
            }
            $self->{other} = $self->get_other( $name, $self->event );
            $self->{other_id} = $self->{other}->dbID if ( $self->{other} );
        }
    }
    return $self->{other_name};
}

=head2 get_other

  Arg [1]   : ReseqTrack::EventComplete
  Arg [2]   : string, name
  Arg [3]   : ReseqTrack::Event
  Function  : fetchs other object for given name and event object
  Returntype: ReseqTrack::Base
  Exceptions: throws if can't use table name to work out which object to fetch
  Example   : 

=cut

sub get_other {
    my ( $self, $name, $event ) = @_;
    if ( $self->table_name eq 'file' || $self->table_name eq 'run_meta_info' ) {
        return $self->get_other_adaptor->fetch_by_name($name);
    }
    elsif ($self->table_name eq 'collection'
        || $self->table_name eq 'input_string' )
    {
        return $self->get_other_adaptor->fetch_by_name_and_type( $name,
            $event->type );
    }
    else {
        throw( "Don't know how to handle " . $self->table_name );
    }
}

=head2 get_other_adaptor

  Arg [1]   : ReseqTrack::EventComplete
  Function  : uses table name to get apropriate object adaptor
  Returntype: ReseqTrack::DBSQL::BaseAdaptor
  Exceptions: throws if not sure which adaptor to fetch for given table name
  Example   : 

=cut

sub get_other_adaptor {
    my ($self) = @_;
    if ( $self->table_name eq "file" ) {
        return $self->adaptor->db->get_FileAdaptor;
    }
    elsif ( $self->table_name eq 'collection' ) {
        return $self->adaptor->db->get_CollectionAdaptor;
    }
    elsif ( $self->table_name eq 'run_meta_info' ) {
        return $self->adaptor->db->get_RunMetaInfoAdaptor;
    }
    elsif ( $self->table_name eq 'input_string' ) {
        return $self->adaptor->db->get_InputStringAdaptor;
    }
    else {
        throw(
            "Don't know what sort of adaptor to get for " . $self->table_name );
    }
}

1;

