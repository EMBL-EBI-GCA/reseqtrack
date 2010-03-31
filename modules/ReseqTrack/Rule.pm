package ReseqTrack::Rule;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Base;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::Base);

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($goal_event, $conditions) = rearrange([qw(GOAL_EVENT CONDITIONS)], @args);
  $self->goal_event($goal_event);
  $self->conditions($conditions);
  return $self;
}


sub goal_event{
  my ($self, $event) = @_;
  if($event){
    throw("Rule::goal_event must be an ReseqTrack::Event object ".
          " not at ".$event) unless($event->isa('ReseqTrack::Event'));
    $self->{'goal_event'} = $event;
  }
  if(!$self->{'goal_event'} && $self->adaptor){
    $self->{'goal_event'} = $self->get_goal_event;
  }
  return $self->{'goal_event'};
}

sub conditions{
  my ($self, $conditions) = @_;
  if($conditions){
    if($conditions->isa('ReseqTrack::Event')){
      push(@{$self->{'conditions'}}, $conditions);
    }elsif(ref($conditions) eq 'ARRAY'){
      throw("Rule::Conditions does not accept array refs of ".$conditions->[0])
          unless($conditions->[0]->isa('ReseqTrack::Event'));  
      push(@{$self->{'conditions'}}, @$conditions);
    }else{
      throw("Must give Rule::conditions an event object or an arrayref not ".$conditions);
    }
  }
  if(!$self->{'conditions'} && $self->adaptor){
    $self->{'conditions'} = $self->get_conditions;
  }
  return $self->{'conditions'};
}

sub get_goal_event{
  my ($self) = @_;
  return $self->adaptor->fetch_goal_event($self->dbID);
}

sub get_conditions{
  my ($self) = @_;
  return $self->adaptor->fetch_conditional_analyses($self->dbID);
}
