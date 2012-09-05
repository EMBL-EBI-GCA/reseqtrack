package ReseqTrack::DBSQL::InputStringAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::InputString;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


#update

sub columns{
  return "input_string.input_string_id, input_string.name, input_string.type";
}

sub table_name{
  return "input_string";
}



sub fetch_by_type{
  my ($self, $type) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name;
  $sql .= " where type = ?";
  my @objects;
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $type);
  $sth->execute();
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@objects, $object);
  }
  $sth->finish;
  return \@objects;
}

sub fetch_incomplete_by_event{
  my ($self, $event) = @_;
  throw($event->name . " table name is not " . $self->table_name) if ($event->table_name ne $self->table_name);
  my $table_name = $self->table_name;
  my $sql = "select ".$self->columns." from $table_name ".
      "left outer join (select event_complete.other_id ".
      "from event_complete where event_complete.event_id = ?) as e ".
      "on $table_name.input_string_id = e.other_id ".
      "where $table_name.type=? and e.other_id is null";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $event->dbID);
  $sth->bind_param(2, $event->type);
  $sth->execute;
  my @input_strings;
  while(my $hashref = $sth->fetchrow_hashref){
    my $input_string = $self->object_from_hashref($hashref);
    push(@input_strings, $input_string);
  }
  
  $sth->finish;
  return \@input_strings;
}

sub fetch_by_name{
  my ($self, $name) = @_;
  my $sql = "select ".$self->columns." from input_string where name = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name);
  $sth->execute;
  my @input_strings;
  $sth->execute;
  while(my ($rowHashref) = $sth->fetchrow_hashref){
    my $input_string = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@input_strings, $input_string);
  }
  $sth->finish;
  return \@input_strings;
}
sub fetch_by_name_and_type{
  my ($self, $name, $type) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name;
  $sql .= " where type = ? and name = ?";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $type);
  $sth->bind_param(2, $name);
  $sth->execute();
  my $rowHashref = $sth->fetchrow_hashref;
  my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;
  return $object;
}


#sub fetch_by_name_and_type{
#  my ($self, $name, $type) = @_;
#  my $sql = "select ".$self->columns." from input_string where name = ? ".
#      "and type = ?";
#  my $sth = $self->prepare($sql);
#  $sth->bind_param(1, $name);
#  $sth->bind_param(2, $type);
#  $sth->execute;
 # my @input_strings;
 # while(my ($rowHashref) = $sth->fetchrow_hashref){
 #   my $input_string = $self->object_from_hashref($rowHashref) if($rowHashref);
 #   push(@input_strings, $input_string);
 # }
 # $sth->finish;
 # return \@input_strings;
#}


sub store{
  my ($self, $input_string) = @_;
  my $sql = "insert ignore into input_string ".
      "(name, type) ".
      "values(?, ?) ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $input_string->name);
  $sth->bind_param(2, $input_string->type);
  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
  if(!$dbID){
    my $new_input_string = $self->fetch_by_name_and_type($input_string->name, 
                                                         $input_string->type);
    $input_string->dbID($new_input_string->dbID);
  }
  $sth->finish();
  $input_string->dbID($dbID) if($dbID);
  $input_string->adaptor($self);
  return $input_string;
}


sub update{
  my ($self, $input_string) = @_;
  throw("Can't update a input_string object without a dbID") unless($input_string->dbID);
  my $sql = "update input_string ".
      "set name = ?, ".
      "type = ?".
      " where input_string_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $input_string->name);
  $sth->bind_param(2, $input_string->type);
  $sth->bind_param(3, $input_string->dbID);
  $sth->execute();
  $sth->finish();
  return $input_string;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create an object from an undefined hashref") if(!$hashref);
  
  my $input_string = ReseqTrack::InputString->new
      (
       -name => $hashref->{name},
       -adaptor => $self,
       -dbID => $hashref->{input_string_id},
       -remote => $hashref->{type},
      );
  return $input_string;
}



1;
