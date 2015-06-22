package ReseqTrack::DBSQL::EventAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Event;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument  qw(rearrange);


@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);




sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


sub columns{

  return  " event.event_id, event.name, event.program,  event.options, event.input_flag, event.farm_options, event.runner_options, event.max_array_size, event.job_slot_limit,  event.output_path, event.type,event.table_name, event.created, event.updated ";
  
}

sub table_name{
  return "event";
}



sub fetch_by_name{

  my ($self, $name) = @_;

  my $sql = "select ".$self->columns." from event where name = ?";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name);
  $sth->execute;
  my ($rowHashref) = $sth->fetchrow_hashref;
  my $event = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;

  return $event;
}
#############


sub store{
  my ($self, $event) = @_;
  my $sql = 
      "INSERT INTO  event " 
      . "(name, program, options,
          input_flag, farm_options, runner_options, max_array_size, job_slot_limit, output_path,
          type, table_name, created, updated) "
      . "values(?, ?, ?, ?, ?, ?, ?, ?, ?,?, ?, now(), now() ) ";
 
  my $sth = $self->prepare($sql);

  $sth->bind_param(1,  $event->name);
  $sth->bind_param(2,  $event->program);
  $sth->bind_param(3,  $event->options);
  $sth->bind_param(4,  $event->input_flag);
  $sth->bind_param(5,  $event->farm_options);
  $sth->bind_param(6,  $event->runner_options);
  $sth->bind_param(7,  $event->max_array_size);
  $sth->bind_param(8,  $event->job_slot_limit);
  $sth->bind_param(9,  $event->output_path);
  $sth->bind_param(10,  $event->type);
  $sth->bind_param(11, $event->table_name);

  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
 
  $sth->finish();
  $event->dbID($dbID);
  $event->adaptor($self);
  $self->store_statistics($event);
  $self->store_history($event);
}
#############


sub update{
  my ($self, $event) = @_;

  # print "Trying to update ". $event->name. "\n";
  
  my $sql = 
      "UPDATE event SET program   = ?, ".
      "options     = ?, ".
      "input_flag = ?, ".
      "farm_options = ?,  ".
      "runner_options = ?,  ".
      "max_array_size = ?,  ".
      "job_slot_limit = ?,  ".
      "output_path = ?, ".
      "type  = ?, ".
      "table_name   = ?,  ".
      "name = ?, ".
      "updated = now() ". 
      "WHERE event_id  = ?;";
  
  
  my $sth = $self->prepare($sql);
  
  
  $sth->bind_param(1,  $event->program);
  $sth->bind_param(2,  $event->options);
  $sth->bind_param(3,  $event->input_flag);
  $sth->bind_param(4,  $event->farm_options);
  $sth->bind_param(5,  $event->runner_options);
  $sth->bind_param(6,  $event->max_array_size);
  $sth->bind_param(7,  $event->job_slot_limit);
  $sth->bind_param(8,  $event->output_path);
  $sth->bind_param(9,  $event->type);
  $sth->bind_param(10,  $event->table_name);
  $sth->bind_param(11, $event->name);
  $sth->bind_param(12, $event->dbID);
 
  $sth->execute();
  $sth->finish();
  $self->store_statistics($event);
  $self->store_history($event);
  return;
}
############


sub object_from_hashref{
    my ($self, $hashref) = @_;

    throw("Can't create an object from an undefined hashref") if(!$hashref);  

    my $event = ReseqTrack::Event->new
	(
         -dbID => $hashref->{event_id},
         -adaptor => $self,
	 -name            =>$hashref->{name},
	 -program         =>$hashref->{program},
	 -options         =>$hashref->{options},
	 -input_flag      =>$hashref->{input_flag},
	 -farm_options    =>$hashref->{farm_options},
	 -runner_options  =>$hashref->{runner_options},
	 -max_array_size  =>$hashref->{max_array_size},
	 -job_slot_limit  =>$hashref->{job_slot_limit},
	 -output_path      =>$hashref->{output_path},
	 -type                 =>$hashref->{type},
	 -table_name     =>$hashref->{table_name},
	 -created           =>$hashref->{created},
	 -updated          =>$hashref->{updated},
     
    );
    return $event; 
}


1;
