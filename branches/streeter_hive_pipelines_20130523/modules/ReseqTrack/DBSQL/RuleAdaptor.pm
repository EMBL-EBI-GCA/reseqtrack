package ReseqTrack::DBSQL::WorkflowAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Workflow;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);


sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


sub columns{
  return "workflow_goal.workflow_id";
}

sub table_name{
  return "workflow_goal";
}


sub fetch_by_goal_event_id{
  my ($self, $goal_event_id) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name." where ".$self->where.
      " and workflow_goal.goal_event_id = ?";
   my $sth = $self->prepare($sql);
   $sth->bind_param(1, $goal_event_id);
   $sth->execute;
   my $rowHashref = $sth->fetchrow_hashref;
   my $object = $self->object_from_hashref($rowHashref);
   $sth->finish;
   return $object;
}

sub fetch_by_goal_event{
  my ($self, $event) = @_;
  return $self->fetch_by_goal_event_id($event->dbID);
}

sub fetch_by_conditional_event_id{
  my ($self, $event_id) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name." where ".$self->where.
      " and workflow_conditions.conditional_event_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $event_id);
  $sth->execute;
  my @objects;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref);
    push(@objects, $object);
  }
  $sth->finish;
  return \@objects;
}

sub fetch_by_conditional_event{
  my ($self, $event) = @_;
  return $self->fetch_by_conditional_event_id($event->dbID);
}

sub fetch_goal_event{
  my ($self, $dbID) = @_;
  my $sql = "select workflow_goal.goal_event_id from workflow_goal where workflow_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $dbID);
  $sth->execute;
  my $rowHashref = $sth->fetchrow_hashref;
  my $aa = $self->db->get_AnalysisAdaptor;
  my $event = $aa->fetch_by_dbID($rowHashref->{goal_event_id});
  $sth->finish;
  return $event;
}

sub fetch_conditional_events{
  my ($self, $dbID) = @_;
  my $sql = "select workflow_conditions.conditional_event_id from workflow_conditions where ".
      "workflow_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $dbID);
  $sth->execute;
  my @events;
  my $aa = $self->db->get_AnalysisAdaptor;
  while (my $rowHashref = $sth->fetchrow_hashref){
    push(@events, $aa->fetch_by_dbID($rowHashref->{conditional_event_id}));
  }
  return(\@events);
}

sub store{
  my ($self, $workflow) = @_;
  throw("Can't store ".$workflow." must be ReseqTrack::Workflow") 
      unless($workflow->isa("ReseqTrack::Workflow"));
  my $goal_sql = "insert ignore into workflow_goal (goal_event_id) values(?)";
  my $goal_sth = $self->prepare($goal_sql);
  $goal_sth->bind_param(1, $workflow->goal_event->dbID);
  my $rows_inserted = $goal_sth->execute;
  my $dbID = $goal_sth->{'mysql_insertid'};
  $workflow->dbID($dbID);
  $workflow->adaptor($self);
  if(!$workflow->dbID){
    throw("Failed to insert ".$workflow." into the workflow_goal table");
  }
  $goal_sth->finish;
  $self->add_condition($workflow->conditions);
}

sub add_condition{
  my ($self, $workflow, $events) = @_;
  my $sql = "insert ignore into workflow_conditions values(?. ?)";
  my $condition_sth = $self->prepare($sql);
  foreach my $event(@$events){
    $condition_sth->bind_param(1, $workflow->dbID);
    $condition_sth->bind_param(2, $event->dbID);
    $condition_sth->execute;
  }
  $condition_sth->finish;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  my $workflow = ReseqTrack::Workflow>new(
    -dbID => $hashref->{workflow_id},
    -adaptor => $self,
      );
  return $workflow;
}
