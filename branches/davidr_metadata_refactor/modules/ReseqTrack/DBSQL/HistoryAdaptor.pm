package ReseqTrack::DBSQL::HistoryAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::History;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub columns{
  return "history.history_id, history.other_id, history.table_name, ".
      "history.comment, history.time";
}

sub table_name{
  return "history";
}

sub fetch_by_table_name{
  my ($self, $table_name) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      "where table_name = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $table_name);
  $sth->execute;
  my @history_objects;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@history_objects, $object) if($object);
  }
  $sth->finish;
  return \@history_objects;
}

sub fetch_by_other_id_and_table_name{
  my ($self, $other_id, $table_name)  = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where table_name = ? and other_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $table_name);
  $sth->bind_param(2, $other_id);
  $sth->execute;
  my @history_objects;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@history_objects, $object) if($object);
  }
  $sth->finish;
  return \@history_objects;
}

sub fetch_by_comment{
  my ($self, $comment) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      "where comment = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $comment);
  $sth->execute;
  my @history_objects;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@history_objects, $object) if($object);
  }
  $sth->finish;
  return \@history_objects;
}

sub store{
  my ($self, $history) = @_;
  return $history if ($history->dbID && 
                      $history->adaptor &&
                      $history->adaptor->dbc->dbname eq $self->dbc->dbname);
  my $sql = "insert into history (other_id, table_name, comment, time) ".
      "values(?, ?, ?, now())";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $history->other_id);
  $sth->bind_param(2, $history->table_name);
  $sth->bind_param(3, $history->comment);
  $sth->execute;
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish;
  $history->dbID($dbID);
  $history->adaptor($self);
  return $history;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create a History object from an empty hashref") unless($hashref);
  my $history = ReseqTrack::History->new(
    -dbID => $hashref->{history_id},
    -adaptor => $self,
    -other_id => $hashref->{other_id},
    -table_name => $hashref->{table_name},
    -comment => $hashref->{comment},
    -time => $hashref->{time},
      );
  return $history;
}

1;
