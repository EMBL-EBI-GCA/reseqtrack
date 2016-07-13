package ReseqTrack::DBSQL::AttributeAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Attribute;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub columns{
  return "attribute.attribute_id, attribute.other_id, attribute.table_name, ".
      "attribute.attribute_name, attribute.attribute_value, attribute.attribute_units";
}

sub table_name{
  return "attribute";
}

sub fetch_by_table_name{
  my ($self, $table_name) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      "where table_name = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $table_name);
  $sth->execute;
  my @attributes_objects;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@attributes_objects, $object) if($object);
  }
  $sth->finish;
  return \@attributes_objects;
}

sub fetch_by_other_id_and_table_name{
  my ($self, $other_id, $table_name)  = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where table_name = ? and other_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $table_name);
  $sth->bind_param(2, $other_id);
  $sth->execute;
  my @attributes_objects;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@attributes_objects, $object) if($object);
  }
  $sth->finish;
  return \@attributes_objects;
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
  my ($self, $attribute, $update) = @_;
  
  throw("Can't store attribute without an other_id") unless($attribute->other_id);
  throw("Can't store attribute without a table_nme") unless($attribute->table_name);
  
  my $exists = $self->fetch_by_other_id_and_table_name_and_attribute_name($attribute->other_id, $attribute->table_name, $attribute->attribute_name);
  
  if($exists){
  	$attribute->dbID($exists->dbID);
  	$attribute->adaptor($self);
  	
    if($update &&
      ($exists->attribute_value ne $attribute->attribute_value || $exists->attribute_units // '' ne $attribute->attribute_units // '')){
      $attribute->dbID($exists->dbID);
#      print "existing stats id is " . $exists->dbID . "\n";
#			print "reassigned stats id is: " . $attribute->dbID . "\n";
#			print "new attribute value is " . $attribute->attribute_value . "\n";
			return $self->update($attribute);      
    }
    return $attribute;
  }
  
  
  my $sql = "insert into attribute (other_id, table_name, attribute_name, ".
      "attribute_value, attribute_units) values(?, ?, ?, ?, ?)";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $attribute->other_id);
  $sth->bind_param(2, $attribute->table_name);
  $sth->bind_param(3, $attribute->attribute_name);
  $sth->bind_param(4, $attribute->attribute_value);
  $sth->bind_param(5, $attribute->attribute_units);
  $sth->execute;
  
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish;
  $attribute->dbID($dbID);
  $attribute->adaptor($self);
  return $attribute;
}


sub update{
  my ($self, $attribute) = @_;

  my $sql = "update attribute ".
      "set table_name = ? ".
      ", other_id = ? ".
      ", attribute_name = ? ".
      ", attribute_value = ? ".
      ", attribute_units = ? ".
      "where attribute_id = ? ";
  my $sth = $self->prepare($sql);

  $sth->bind_param(1, $attribute->table_name);
  $sth->bind_param(2, $attribute->other_id);
  $sth->bind_param(3, $attribute->attribute_name);
  $sth->bind_param(4, $attribute->attribute_value);
  $sth->bind_param(5, $attribute->attribute_units);
  $sth->bind_param(6, $attribute->dbID);  
  $sth->execute; 
  $sth->finish; 
  
  return $attribute; 
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create a History object from an empty hashref") unless($hashref);
  #print "test object_from_hashref, stats id is " . $hashref->{attribute_id} . "\n";
  my $attribute = ReseqTrack::Attribute->new(
    -dbID =>  $hashref->{attribute_id},
    -adaptor => $self,
    -other_id => $hashref->{other_id},
    -table_name => $hashref->{table_name},
    -attribute_name => $hashref->{attribute_name},
    -attribute_value => $hashref->{attribute_value},
    -attribute_units => $hashref->{attribute_units},
      );
  return $attribute;
}

sub delete {
	my ($self, $attribute) = @_;
	my $sql = "delete from ".$self->table_name()." where attribute_id = ?";
	my $sth = $self->prepare($sql);
	$sth->bind_param(1, $attribute->dbID);
	$sth->execute();
	$sth->finish();
}

1;
