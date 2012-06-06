package ReseqTrack::DBSQL::StatisticsAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Statistic;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub columns{
  return "statistics.statistics_id, statistics.other_id, statistics.table_name, ".
      "statistics.attribute_name, statistics.attribute_value";
}

sub table_name{
  return "statistics";
}

sub fetch_by_table_name{
  my ($self, $table_name) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      "where table_name = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $table_name);
  $sth->execute;
  my @statistics_objects;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@statistics_objects, $object) if($object);
  }
  $sth->finish;
  return \@statistics_objects;
}

sub fetch_by_other_id_and_table_name{
  my ($self, $other_id, $table_name)  = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where table_name = ? and other_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $table_name);
  $sth->bind_param(2, $other_id);
  $sth->execute;
  my @statistics_objects;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@statistics_objects, $object) if($object);
  }
  $sth->finish;
  return \@statistics_objects;
}

sub fetch_by_other_id_and_table_name_and_attribute_name{
  my ($self, $other_id, $table_name, $attribute_name)  = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where table_name = ? and other_id = ? and attribute_name = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $table_name);
  $sth->bind_param(2, $other_id);
  $sth->bind_param(3, $attribute_name);
  $sth->execute;
  my $rowHashref = $sth->fetchrow_hashref;
  my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;
  #print "fetched stats obj stats id: " . $object->dbID . "\n";
  return $object;
}

sub store{
  my ($self, $statistics, $update) = @_;
  my $exists = $self->fetch_by_other_id_and_table_name_and_attribute_name
      ($statistics->other_id, $statistics->table_name, $statistics->attribute_name);
  if($exists){
    if($exists->attribute_value eq $statistics->attribute_value){
      $statistics->dbID($exists->dbID);
      $statistics->adaptor($self);
      return $statistics;
    }else{
      if($update){
	$statistics->dbID($exists->dbID); #holly added this line
       	#print "existing stats id is " . $exists->dbID . "\n";
	#print "reassigned stats id is: " . $statistics->dbID . "\n";
	#print "new attribute value is " . $statistics->attribute_value . "\n";
	return $self->update($statistics);
      }
    }
  }
  throw("Can't store statistic without an other id") unless($statistics->other_id);
  my $sql = "insert into statistics (other_id, table_name, attribute_name, ".
      "attribute_value) values(?, ?, ?, ?)";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $statistics->other_id);
  $sth->bind_param(2, $statistics->table_name);
  $sth->bind_param(3, $statistics->attribute_name);
  $sth->bind_param(4, $statistics->attribute_value);
  $sth->execute;
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish;
  $statistics->dbID($dbID);
  $statistics->adaptor($self);
  return $statistics;
}


sub update{
  my ($self, $statistics) = @_;

  my $sql = "update statistics ".
      "set table_name = ? ".
      ", other_id = ? ".
      ", attribute_name = ? ".
      ", attribute_value = ? ".
      "where statistics_id = ? ";
  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $statistics->table_name);
  $sth->bind_param(2, $statistics->other_id);
  $sth->bind_param(3, $statistics->attribute_name);
  $sth->bind_param(4, $statistics->attribute_value);
  $sth->bind_param(5, $statistics->dbID);  
  $sth->execute; 
  $sth->finish; 
  return $statistics; 
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create a History object from an empty hashref") unless($hashref);
  #print "test object_from_hashref, stats id is " . $hashref->{statistics_id} . "\n";
  my $statistics = ReseqTrack::Statistic->new(
    #-db_id => $hashref->{statistics_id},
    -dbID =>  $hashref->{statistics_id},
    -adaptor => $self,
    -other_id => $hashref->{other_id},
    -table_name => $hashref->{table_name},
    -attribute_name => $hashref->{attribute_name},
    -attribute_value => $hashref->{attribute_value},
      );
  return $statistics;
}

1;
