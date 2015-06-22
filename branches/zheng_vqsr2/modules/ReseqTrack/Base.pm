
=pod

=head1 NAME

ReseqTrack::Base

=head1 SYNOPSIS

This is a base class for database objects to provide the standard adaptor and dbID
methods

This object forms the base of the inheritance tree

=head1 Example


=cut

package ReseqTrack::Base;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

sub new {
    my $caller = shift;
    my $class = ref($caller) || $caller;

    my ( $adaptor, $dbID ) = rearrange( [ 'ADAPTOR', 'dbID' ], @_ );

    if ($adaptor) {
        if (   !ref($adaptor)
            || !$adaptor->isa('ReseqTrack::DBSQL::BaseAdaptor') )
        {
            throw('-ADAPTOR argument must be a ReseqTrack::DBSQL::BaseAdaptor');
        }
    }
    return bless( { 'dbID' => $dbID, 'adaptor' => $adaptor }, $class );
}

=head2 dbID

  Arg [1]    : int $dbID
  Example    : none
  Description: getter/setter for the database internal id
  Returntype : int
  Exceptions : none
  Caller     : general, set from adaptor on store
  Status     : Stable

=cut

sub dbID {
    my $self = shift;
    $self->{'dbID'} = shift if (@_);
    return $self->{'dbID'};
}

=head2 adaptor

  Arg [1]    : Bio::EnsEMBL::DBSQL::BaseAdaptor $adaptor
  Example    : none
  Description: get/set for this objects Adaptor
  Returntype : Bio::EnsEMBL::DBSQL::BaseAdaptor
  Exceptions : none
  Caller     : general, set from adaptor on store
  Status     : Stable

=cut

sub adaptor {
    my $self = shift;

    if (@_) {
        my $ad = shift;
        if ( $ad
            && ( !ref($ad) || !$ad->isa('ReseqTrack::DBSQL::BaseAdaptor') ) )
        {
            throw('Adaptor argument must be a ReseqTrack::DBSQL::BaseAdaptor');
        }
        $self->{'adaptor'} = $ad;
    }

    return $self->{'adaptor'};
}

=head2 object_table_name

  Arg [1]   : ReseqTrack::Base
  Function  : this is a method which should be implemented by all
  child objects of Base, the method should return a string of the
  table name which holds that object
  Returntype: string
  Exceptions: throws as this means the object isn't implemented
  Example   : 

=cut

sub object_table_name {
    my ($self) = @_;
    throw(  "ReseqTrack::Base: need to implement an object table name"
          . " method in "
          . $self );
}

1;
