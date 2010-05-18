package ReseqTrack::DBSQL::CollectionAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Collection;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use File::Basename;

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub table_name{
  return "collection";
}

sub columns{
  return "collection.collection_id, collection.name, collection.type, collection.table_name";
}

sub fetch_by_name{
  my ($self, $name) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where name = ?";
  my @collections;
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name);
  $sth->execute;
  while(my $hashref = $sth->fetchrow_hashref){
    my $collection = $self->object_from_hashref($hashref) if($hashref);
    push(@collections, $collection) if($collection);
  }
  $sth->finish;
  return \@collections;
}

sub fetch_by_type{
  my ($self, $type) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where type = ?";
  my @collections;
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $type);
  $sth->execute;
  while(my $hashref = $sth->fetchrow_hashref){
    my $collection = $self->object_from_hashref($hashref) if($hashref);
    push(@collections, $collection) if($collection);
  }
  
  $sth->finish;
  return \@collections;
}

sub fetch_by_name_and_type{
  my ($self, $name, $type) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where name = ? and type = ?";
  #print $sql."\n";
  my @collections;
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name);
  $sth->bind_param(2, $type);
  $sth->execute;
  my $hashref = $sth->fetchrow_hashref;
  my $collection = $self->object_from_hashref($hashref) if($hashref);
  $sth->finish;
  return $collection;
}

sub fetch_by_name_and_table_name{
  my ($self, $name, $table_name) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where name = ? and table_name = ?";
  #print $sql."\n";
  my @collections;
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name);
  $sth->bind_param(2, $table_name);
  $sth->execute;
  while(my $hashref = $sth->fetchrow_hashref){
    my $collection = $self->object_from_hashref($hashref) if($hashref);
    push(@collections, $collection);
  }
  $sth->finish;
  return \@collections;
}

sub fetch_by_other_id_and_type{
  my ($self, $other_id, $type) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.", ".
      "collection_group ".
      "where collection.collection_id = collection_group.collection_id ".
      " and collection.type = ? and collection_group.other_id = ?";
  my @collections;
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $type);
  $sth->bind_param(2, $other_id);
  $sth->execute;
  while(my $hashref = $sth->fetchrow_hashref){
    my $collection = $self->object_from_hashref($hashref) if($hashref);
    push(@collections, $collection);
  }
  $sth->finish;
  return \@collections;
}

sub fetch_other_ids{
  my ($self, $collection) = @_;
  my $sql = "select collection_group.other_id ".
      "from collection_group ".
      "where collection_group.collection_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $collection->dbID);
  $sth->execute;
  my @other_ids;
  while(my $hashref = $sth->fetchrow_hashref){
    push(@other_ids, $hashref->{other_id});
  }
  $sth->finish;
  #warning($collection->dbID." ".$collection->name." has no other ids")
  #    if(!@other_ids || @other_ids == 0);
  return \@other_ids;
}

sub store{
  my ($self, $collection, $update) = @_;
  my $others = $collection->others;
  throw("Can't store an empty collection, it has no other objects attached to it")
      unless($others && @$others >= 1);
  #Must check for the existence of the collection with that name and type first
  #print "Checking in ".$self->dbc->dbname."\n";
  my $exists = $self->fetch_by_name_and_type($collection->name, $collection->type);
  if($exists){
    warning($exists->name." ".$exists->type." is already in the datbase skipping");
    $collection->dbID($exists->dbID);
    $collection->adaptor($self);
    $self->store_statistics($collection, $update);
    $self->store_others($collection);
    return $collection;
  }
  throw("Can't store a collection ".$collection." without a table name ".
        $collection->table_name) unless($collection->table_name);
  my $sql = "insert into collection(name, type, table_name) values(?, ?, ?)";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $collection->name);
  $sth->bind_param(2, $collection->type);
  $sth->bind_param(3, $collection->table_name);
  $sth->execute;
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish;
  $collection->dbID($dbID);
  $collection->adaptor($self);
  $self->store_statistics($collection, $update);
  $self->store_others($collection);
  return $collection;
}


sub store_others{
  my ($self, $collection) = @_;
  my $others = $collection->others;
  my $sql = "insert ignore into collection_group(collection_id, other_id) values(?, ?)";
  my $sth = $self->prepare($sql);
  my $oa = $collection->get_other_adaptor;
  my $existing_others = $self->fetch_other_ids($collection);
  my %existing_hash;
  foreach my $other_id(@$existing_others){
    $existing_hash{$other_id} = 1;
  }
  foreach my $other(@$others){
    unless($other->dbID){
      $other = $oa->store($other);
      throw($oa.":store has failed to generate a dbID for the object this ".
            "means we can't associate it with the collection ") unless($other->dbID);
    }
    next if($existing_hash{$other->dbID});
    $sth->bind_param(1, $collection->dbID);
    $sth->bind_param(2, $other->dbID);
    $sth->execute;
  }
  $sth->finish;
}



sub update_type{
  my ($self, $collection) = @_;
  if(!$collection->dbID){
    throw("Can't update a collection object ".$collection->name." ".
          $collection->type." without a dbID");
  }
  if(!$collection->history || @{$collection->history} == 0){
    throw("Can't update a collection object without a history object");
  }
  my $sql = "update collection ".
      "set type = ?".
      "where collection_id = ?";
  my $sth = $self->prepare($sql);
  $sth->execute($collection->type, $collection->dbID);
  $self->store_history($collection);
  return $collection;
}
sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create ReseqTrack::Collection from an empty hashref")
      if(!$hashref || keys(%$hashref) == 0);
  # print $hashref->{collection_id}." ".$hashref->{table_name}."\n";
  throw("Can't create a ReseqTrack::Collection based on ".$hashref->{collection_id}.
        " without a table name") unless($hashref->{table_name});
  my $collection = ReseqTrack::Collection->new(
    -dbID => $hashref->{collection_id},
    -adaptor => $self,
    -name => $hashref->{name},
    -type => $hashref->{type},
    -table_name => $hashref->{table_name}
      );
  return $collection;
}




1;
