
=pod

=head1 NAME

ReseqTrack::ArchiveAction

=head1 SYNOPSIS

This is a container object for the archive_action table which provides human memorable
names for the different archive action ids

=cut

package ReseqTrack::ArchiveAction;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

=head2 new

  Arg [1]   : ReseqTrack::ArchiveAction
  Arg [2]   : string, name of the action
  Function  : create ReseqTrack::ArchiveAction object
  Returntype: ReseqTrack::ArchiveAction
  Exceptions: throws if dbID or action name aren't recognised
  Example   : 

=cut

sub new {
   my ( $class, @args ) = @_;
   my $self = $class->SUPER::new(@args);
   my ($action) = rearrange( ['ACTION'], @args );
   $self->action($action);
   throw( "Can't have a dbID " . $self->dbID . " which isn't 1, 2, 3, 4, 5 or 6" )
     unless ( $self->dbID == 1
       || $self->dbID == 2
       || $self->dbID == 3
       || $self->dbID == 4
       || $self->dbID == 5
       || $self->dbID == 6 );
   throw(  "Can't have an action "
         . $self->action
         . " which isn't archive, dearchive, "
         . "replace or move_within_volume" )
     unless ( $self->action eq "archive"
       || $self->action eq "dearchive"
       || $self->action eq "replace"
       || $self->action eq "move_within_volume"
       || $self->action eq "delete_permanently"
       || $self->action eq "copy_to_staging");
   return $self;
}


=head2 action

  Arg [1]   : ReseqTrack::ArchiveAction
  Arg [2]   : string, action name
  Function  : accessor method for action name
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub action {
    my ( $self, $action ) = @_;
    if ($action) {
        $self->{'action'} = $action;
    }
    return $self->{'action'};
}

=head2 object_table_name

  Arg [1]   : ReseqTrack::ArchiveAction
  Function  : returning table name
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub object_table_name {
    my ($self) = @_;
    return "archive_action";
}

1;
