
=pod

=head1 NAME

ReseqTrack::Collection

=head1 SYNOPSIS

This object is designed to represent groups of other objects from the run_meta_info,
 file or alignment_meta_info tables

=head1 Example

my $collection =  ReseqTrack::Collection->new(
  -name => 'ERR000270',
  -others => \@file_objects,
  -type => $file_objects[0]->type
  -table_name => 'file',
    );

=cut

package ReseqTrack::Collection;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::HasHistory;

@ISA = qw(ReseqTrack::HasHistory);

=head2 new

  Arg [1]   : ReseqTrack::Collection
  Arg [2]   : string, name
  Arg [3]   : string, type
  Arg [4]   : arrayref, of other objects
  Arg [5]   : arrayref, of other object dbIDs
  Arg [6]   : string, table name
  Function  : create a ReseqTrack::Collection object
  Returntype: ReseqTrack::Collection

  Exceptions: throws if no name, type or table name is defined
  throws if the table name is not one of the accepted types
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my ( $name, $type, $others, $other_ids, $table_name ) =
      rearrange( [qw(NAME TYPE OTHERS OTHER_IDS TABLE_NAME)], @args );
    ####error checking

    throw("Can't create a ReseqTrack::Collection object without a name")
      unless ($name);
    throw("Can't create a ReseqTrack::Collection object without a type")
      unless ($type);
    throw("Can't create a ReseqTrack::Collection object wihtout a table_name")
      unless ($table_name);
    if ($table_name) {
        throw(
"ReqseqTrack::Collection::table_name must be file, alignment_meta_info "
              . "run_meta_info or collection not "
              . $table_name )
          unless ( $table_name eq 'file'
            || $table_name eq 'alignment_meta_info'
            || $table_name eq 'collection'
            || $table_name eq 'run_meta_info' );
    }
    #########

    $self->name($name);
    $self->type($type);
    $self->other_ids($other_ids);
    $self->others($others);
    $self->table_name($table_name);
    return $self;
}

=head2 Accessor methods

  Arg [1]   : ReseqTrack::Collection
  Arg [2]   : int/string
  Function  : accessor method for variables
  Returntype: int/string
  Exceptions: 
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
    if ($arg) {
        $self->{type} = $arg;
    }
    return $self->{type};
}

=head2 other_ids

  Arg [1]   : ReseqTrack::Collection
  Arg [2]   : arrayref of other object dbIDs
  Function  : store the array of other object dbIDs and fetch list of other 
  ids from the database is there is no list and the adpator is defined
  Returntype: arrayref
  Exceptions: 
  Example   : 

=cut

sub other_ids {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{other_ids} = $arg;
    }
    if ( !$self->{other_ids} ) {
        if ( $self->adaptor ) {
            my $ids = $self->adaptor->fetch_other_ids($self);
            $self->{other_ids} = $ids;
        }
    }
    return $self->{other_ids};
}

=head2 table_name

  Arg [1]   : ReseqTrack::Collection
  Arg [2]   : string, table name
  Function  : store given table name, work out table name based on other objects
  if none defined
  Returntype: 
  Exceptions: 
  Example   : 

=cut

sub table_name {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{table_name} = $arg;
    }
    if ( !$self->{table_name} ) {
        $self->calculate_table_name;
    }
    return $self->{table_name};
}

=head2 calculate_table_name

  Arg [1]   : ReseqTrack::Collection
  Function  : work out table name based on attached other objects
  Returntype: string
  Exceptions: n/a
  Example   : 

=cut

sub calculate_table_name {
    my ($self) = @_;
    if ( $self->{others} ) {
        my $first = $self->{others}->[0];
        $self->{table_name} = $first->object_table_name;
    }
    return $self->{table_name};
}

=head2 get_objects

  Arg [1]   : ReseqTrack::Collection
  Function  : fetch other objects from database
  Returntype: arrayref of other objects
  Exceptions: throws if can't find objects using the specific dbIDs
  Example   : 

=cut

sub get_objects {
    my ($self) = @_;
    $self->sanity_checked_others(0);
    my $adaptor = $self->get_other_adaptor;
    my @objects;
    foreach my $id ( @{ $self->other_ids } ) {
        my $object = $adaptor->fetch_by_dbID($id);
        throw( "Failed to fetch object using " . $adaptor . " and " . $id )
          unless ($object );
	push( @objects, $object );
    }
    $self->{others} = \@objects;
}

=head2 get_other_adaptor

  Arg [1]   : ReseqTrack::Collection
  Function  : gets an object adaptor given the table name
  Returntype: ReseqTrack::DBSQL::BaseAdaptor
  Exceptions: throws if no table name defined if if table name is not recognised
  Example   : 

=cut

sub get_other_adaptor {
    my ($self) = @_;
    if ( !$self->table_name ) {
        throw(
"ReseqTrack::Collection:Don't know what adaptor to get as undefined table name"
        );
    }
    if ( $self->table_name eq "file" ) {
        return $self->adaptor->db->get_FileAdaptor;
    }
    elsif ( $self->table_name eq 'run_meta_info' ) {
        return $self->adaptor->db->get_RunMetaInfoAdaptor;
    }
    elsif ( $self->table_name eq 'alignment_meta_info' ) {
        return $self->adaptor->db->get_AlignmentMetaInfoAdaptor;
    }
    elsif ( $self->table_name eq 'collection' ) {
        return $self->adaptor->db->get_CollectionAdaptor;
    }
    else {
        throw(
            "Don't know what sort of adaptor to get for " . $self->table_name );
    }
}

=head2 others

  Arg [1]   : ReseqTrack::Collection
  Arg [2]   : arrayref or ReseqTrack::Base object
  Function  : holder for arrayref of objects, if no array defined will use other information to fetch objects
  Returntype: arrayref of ReseqTrack::Base objects
  Exceptions: throws if array if of mixed type, ie one is a ReseqTrack::File and another if ReseqTrack::RunMetaInfo
  Example   : 

=cut

sub others {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->sanity_checked_others(0);
        if ( ref($arg) eq 'ARRAY' ) {
            push( @{ $self->{others} }, @$arg );
        }
        else {
            push( @{ $self->{others} }, $arg );
        }
    }
    if ( !$self->{others} ) {
        if ( $self->adaptor && $self->other_ids && $self->table_name ) {
            $self->get_objects;
        }
    }

    unless ( $self->sanity_checked_others ) {
        my $others = $self->{others};
        if ( $others && @$others >= 2 ) {
            my $first = ref( $others->[0] );
            foreach my $object (@$others) {
                throw(  "ReseqTrack::Collection " 
                      . $object
                      . " needs to be "
                      . $first . " not "
                      . "what it is " )
                  unless ( $first eq ref($object) );
            }
        }
        $self->sanity_checked_others(1);
    }
    return $self->{others};
}

=head2 sanity_checked_others

  Arg [1]   : ReseqTrack::Collection
  Arg [2]   : 0/1 binary flag
  Function  : indicates if the array of other objects has been sanity checked
  Returntype: 0/1 binary flag
  Exceptions: 
  Example   : 

=cut

sub sanity_checked_others {
    my ( $self, $arg ) = @_;
    if ( defined($arg) ) {
        $self->{sanity_checked_others} = $arg;
    }
    return $self->{sanity_checked_others};
}

=head2 object_table_name

  Arg [1]   : ReseqTrack::Archive
  Function  : return table name for object, archive
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub object_table_name {
    my ($self) = @_;
    return 'collection';
}

1;
