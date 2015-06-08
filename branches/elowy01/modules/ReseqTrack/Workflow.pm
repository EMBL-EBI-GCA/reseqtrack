
=pod

=head1 NAME

ReseqTrack::Workflow;

=head1 SYNOPSIS

An container object to hold the data workflow_goal and workflow_condition tables.

These tables describe dependencies between different event objects to form
a pipeline

=head1 EXAMPLE

 $workflow = ReseqTrack::Workflow->new(
        -goal_event => $event,
        -conditions => \@conditional_events,
          );
 

=cut

package ReseqTrack::Workflow;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Base;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::Base);

=head2 new

  Arg [GOAL_EVENT] : ReseqTrack::Event the goal event object which must be executed 
  when the workflow can be run
  Arg [CONDITIONS] : arrayref of ReseqTrack::Workflow objects which need to be 
  complete before the goal event can be run
  Function  : create ReseqTrack::Workflow object
  Returntype: ReseqTrack::Workflow
  Exceptions: 
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my ( $goal_event, $conditions ) =
      rearrange( [qw(GOAL_EVENT CONDITIONS)], @args );
    $self->goal_event($goal_event);
    $self->conditions($conditions);
    return $self;
}

=head2 goal_event

  Arg [1]   : ReseqTrack::Event 
  Function  : hold and return event object
  Returntype: ReseqTrack::Event
  Exceptions: throws if passed something other than ReseqTrack::Event
  Example   : 

=cut

sub goal_event {
    my ( $self, $event ) = @_;
    if ($event) {
        throw(  "Workflow::goal_event must be an ReseqTrack::Event object "
              . " not at "
              . $event )
          unless ( $event->isa('ReseqTrack::Event') );
        $self->{'goal_event'} = $event;
    }
    if ( !$self->{'goal_event'} && $self->adaptor ) {
        $self->{'goal_event'} = $self->get_goal_event;
    }
    return $self->{'goal_event'};
}

=head2 conditions

  Arg [1]   : Arrayref of ReseqTrack::Event objects
  Function  : hold and return arrayref of conditional event objects
  Returntype: arrayref of ReseqTrack::Event objects
  Exceptions: throws if not given ReseqTrack::Event objects
  Example   : 

=cut

sub conditions {
    my ( $self, $conditions ) = @_;
    if ($conditions) {
        if ( $conditions->isa('ReseqTrack::Event') ) {
            push( @{ $self->{'conditions'} }, $conditions );
        }
        elsif ( ref($conditions) eq 'ARRAY' ) {
            throw( "Workflow::Conditions does not accept array refs of "
                  . $conditions->[0] )
              unless ( $conditions->[0]->isa('ReseqTrack::Event') );
            push( @{ $self->{'conditions'} }, @$conditions );
        }
        else {
            throw(
"Must give Workflow::conditions an event object or an arrayref not "
                  . $conditions );
        }
    }
    if ( !$self->{'conditions'} && $self->adaptor ) {
        $self->{'conditions'} = $self->get_conditions;
    }
    return $self->{'conditions'};
}

=head2 get_goal_event

  Arg [1]   : n/a
  Function  : fetches goal event object from database 
  Returntype: ReseqTrack::Event
  Exceptions: 
  Example   : 

=cut

sub get_goal_event {
    my ($self) = @_;
    return $self->adaptor->fetch_goal_event( $self->dbID );
}

=head2 get_conditions

  Arg [1]   : n/a
  Function  : get array of conditional events
  Returntype: arrayref of ReseqTrack::Event objects
  Exceptions: 
  Example   : 

=cut

sub get_conditions {
    my ($self) = @_;
    return $self->adaptor->fetch_conditional_events( $self->dbID );
}

1;
