package ReseqTrack::DBSQL::MetaAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;


use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}



sub columns{
  return "meta.meta_id, meta.meta_key, meta.meta_value";
}

sub table_name{
  return "meta";
}


#fetch all
#fetch by key

sub fetch_meta_value_by_meta_key{
  my ($self, $key) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where meta_key = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $key);
  $sth->execute();
  my ($dbID, $db_key, $value) = $sth->fetchrow;
  return $value;
}

#fetch by value
#store

sub store{
  my ($self, $key, $value) = @_;
  my $sql = "insert into meta (meta_key, meta_value) values(?, ?)";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $key);
  $sth->bind_param(2, $value);
  $sth->execute;
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish();
  return $dbID;
}

#update
#remove

sub remove_by_meta_key{
  my ($self, $key) = @_;
  my $sql = "delete from meta where meta_key = ?";
  #print "Deleting ".$key."\n";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $key);
  $sth->execute;
  $sth->finish();
}




1;
