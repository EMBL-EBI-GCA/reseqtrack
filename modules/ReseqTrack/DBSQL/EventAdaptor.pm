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

  return  " event.event_id, event.name, event.program,  event.program_version, event.batch_size, event.options, event.input_flag, event.farm_options, event.output_path, event.type,event.table_name, event.created, event.updated ";
  
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
      . "(name, program, program_version, options,
          input_flag, farm_options, batch_size, output_path,
          type, table_name, created, updated) "
      . "values(?, ?, ?, ?, ?, ?, ?, ?,?, ?, now(), now() ) ";
 
  my $sth = $self->prepare($sql);

  $sth->bind_param(1,  $event->name);
  $sth->bind_param(2,  $event->program);
  $sth->bind_param(3,  $event->program_version);
  $sth->bind_param(4,  $event->options);
  $sth->bind_param(5,  $event->input_flag);
  $sth->bind_param(6,  $event->farm_options);
  $sth->bind_param(7,  $event->batch_size);
  $sth->bind_param(8,  $event->output_path);
  $sth->bind_param(9,  $event->type);
  $sth->bind_param(10, $event->table_name);

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
      "program_version = ?, ".
      "options     = ?, ".
      "input_flag = ?, ".
      "farm_options = ?,  ".
      "batch_size = ?, ". 
      "output_path = ?, ".
      "type  = ?, ".
      "table_name   = ?,  ".
      "name = ?, ".
      "updated = now() ". 
      "WHERE event_id  = ?;";
  
  
  my $sth = $self->prepare($sql);
  
  
  $sth->bind_param(1,  $event->program);
  $sth->bind_param(2,  $event->program_version);
  $sth->bind_param(3,  $event->options);
  $sth->bind_param(4,  $event->input_flag);
  $sth->bind_param(5,  $event->farm_options);
  $sth->bind_param(6,  $event->batch_size);
  $sth->bind_param(7,  $event->output_path);
  $sth->bind_param(8,  $event->type);
  $sth->bind_param(9,  $event->table_name);
  $sth->bind_param(10, $event->name);
  $sth->bind_param(11, $event->dbID);
 
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
	 -program_version =>$hashref->{program_version},
	 -options         =>$hashref->{options},
	 -input_flag      =>$hashref->{input_flag},
	 -farm_options    =>$hashref->{farm_options},
	 -batch_size        =>$hashref->{batch_size},
	 -output_path      =>$hashref->{output_path},
	 -type                 =>$hashref->{type},
	 -table_name     =>$hashref->{table_name},
	 -created           =>$hashref->{created},
	 -updated          =>$hashref->{updated},
     
    );
    return $event; 
}


1;
