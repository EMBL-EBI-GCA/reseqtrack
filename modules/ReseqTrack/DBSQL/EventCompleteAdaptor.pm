package ReseqTrack::DBSQL::EventCompleteAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::EventComplete;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub columns{
  return "event_complete.event_id, event_complete.other_id, ".
      "event_complete.table_name, event_complete.type, event_complete.success, ".
      "event_complete.time";
}

sub table_name{
  return "event_complete";
}


sub fetch_by_dbID{
  throw("ReseqTrack::DBSQL::EventComplete does not have an fetch by dbID method ".
        " as there are no internal indentifiers");
}

sub name_column_name{
  my ($self, $event) = @_;
  if($event->table_name eq 'run_meta_info'){
    return "run_meta_info.run_id";
  }else{
    return $event->table_name.".name";
  }
}

sub fetch_by_input_string_and_event{
  my ($self, $input_string, $event) = @_;
  my $name_column_name = $self->name_column_name($event);
  my $sql = "select ".$self->columns." ".
      "from ".$self->table_name.", ".$event->table_name." ".
      "where ".$self->table_name.".event_id = ? ".
      "and ".$self->table_name.".other_id = ".
      $event->table_name.".".$event->table_name."_id ".
      "and ".$name_column_name." = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $event->dbID);
  $sth->bind_param(2, $input_string);
  $sth->execute;
  my $rowhashref = $sth->fetchrow_hashref;
  my $eventcomplete = $self->object_from_hashref($rowhashref) if($rowhashref);
  return $eventcomplete;
}

sub fetch_by_event{
  my ($self, $event) = @_;
  my $sql = "select ".$self->columns." ".
      "from ".$self->table_name.", ".$event->table_name." ".
      "where ".$self->table_name.".event_id = ? ".
      "and ".$self->table_name.".other_id = ".
      $event->table_name.".".$event->table_name."_id ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $event->dbID);
  $sth->execute;
  my @complete;
  while(my $rowhashref = $sth->fetchrow_hashref){
    my $eventcomplete = $self->object_from_hashref($rowhashref) if($rowhashref);
    push(@complete, $eventcomplete);
  }
  return \@complete;
}


sub fetch_by_other_id{
  my ($self, $other_id, $table_name) = @_;
  my $sql = "select ".$self->columns." ".
      "from ".$self->table_name." ".
      "where other_id = ? ".
      "and table_name =  ? ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $other_id);
  $sth->bind_param(2, $table_name);
  $sth->execute;
  my @complete;
  while(my $rowhashref = $sth->fetchrow_hashref){
    my $eventcomplete = $self->object_from_hashref($rowhashref) if($rowhashref);
    push(@complete, $eventcomplete);
  }
  return \@complete;
}

sub fetch_by_other{
  my ($self, $other, $table_name) = @_;
  if(!$table_name){
    if($other->isa("ReseqTrack::File")){
      $table_name = 'file';
    }elsif($other->isa("ReseqTrack::Collection")){
      $table_name = 'collection';
    }else{
      throw("Can't work out other table name from ".$other);
    }
  }
  return $self->fetch_by_other_id($other, $table_name);
}


sub store{
  my ($self, $event_complete) = @_;
  my $sql = "insert into event_complete(event_id, other_id, table_name, ".
      "type, success, time) values(?, ?, ?, ?, ?, now())";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $event_complete->event_id);
  $sth->bind_param(2, $event_complete->other_id);
  $sth->bind_param(3, $event_complete->table_name);
  $sth->bind_param(4, $event_complete->type);
  $sth->bind_param(5, $event_complete->success);
  $sth->execute;
  $sth->finish;
  $event_complete->adaptor($self);
}

sub remove{
  my ($self, $event_complete) = @_;
  my $sql = "delete from event_complete where event_id = ? and other_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $event_complete->event_id);
  $sth->bind_param(2, $event_complete->other_id);
  $sth->execute;
  $sth->finish;
}

sub set_success{
  my ($self, $event_complete) = @_;
  my $sql = "update event_complete ".
      "set success = ? ".
      "where event_complete.name = ? ";
  my $sth = $self->db->prepare($sql);
  $sth->bind_param(1, $event_complete->name);
  $sth->execute;
  $sth->finish;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create an object from an undefined hashref") if(!$hashref);  
  my $event_complete = ReseqTrack::EventComplete->new
      (
       -adaptor => $self,
       -table_name => $hashref->{table_name},
       -event_id => $hashref->{event_id},
       -other_id => $hashref->{other_id},
       -type => $hashref->{type},
       -success => $hashref->{success},
       -time => $hashref->{time},
      );
  return $event_complete;
}
