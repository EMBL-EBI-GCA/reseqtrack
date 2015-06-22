package ReseqTrack::DBSQL::ArchiveLocationAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::ArchiveLocation;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


sub columns{
  my ($self) = @_;
  return "archive_location.archive_location_id, archive_location.location, archive_location.location_name";
}

sub table_name{
  return "archive_location";
}

sub fetch_by_archive_location{
  my ($self, $location) = @_;
  my $sql = "select ".$self->columns." from archive_location where location = ?";
  my $sth = $self->db->dbc->prepare($sql);
  $sth->bind_param(1, $location);
  $sth->execute();
  my $rowHashref = $sth->fetchrow_hashref;
  my $archive_location = $self->object_from_hashref($rowHashref) if($rowHashref);
  if(!$archive_location){
    if($location =~ /\/$/){
      $location =~ s/\///;
    }else{
      $location .= "/";
    }
    $sth->bind_param(1, $location);
    $sth->execute();
    my $rowHashref = $sth->fetchrow_hashref;
    $archive_location = $self->object_from_hashref($rowHashref) if($rowHashref);
  }
  $sth->finish;
  return $archive_location;
}

sub fetch_by_archive_location_name{
  my ($self, $location_name) = @_;
  my $sql = "select ".$self->columns." from archive_location where location_name = ?";
  my $sth = $self->db->dbc->prepare($sql);
  $sth->bind_param(1, $location_name);
  $sth->execute();
  my $rowHashref = $sth->fetchrow_hashref;
  my $archive_location = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;
  return $archive_location;
}
sub fetch_by_dbID{
  my ($self, $dbID) = @_;
  my $sql = "select ".$self->columns." from archive_location ".
      "where archive_location_id = ?";
  my $sth = $self->db->dbc->prepare($sql);
  $sth->bind_param(1, $dbID);
  $sth->execute();
  my $rowHashref = $sth->fetchrow_hashref;
  my $location = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;
  return $location;
}

sub store{
  my ($self, $archivelocation) = @_;
  my $sql = "insert into archive_location (archive_location_id, location) ".
      " values (?, ?)";
  my $sth = $self->db->dbc->prepare($sql);
  $sth->bind_param(1, $archivelocation->dbID);
  $sth->bind_param(2, $archivelocation->location);
  $sth->execute();
  $sth->finish;
  $archivelocation->adaptor($self);
  return $archivelocation;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create an object for an undefined hashref") if(!$hashref);
    my $location = ReseqTrack::ArchiveLocation->new
      (
       -dbID => $hashref->{archive_location_id},
       -location => $hashref->{location},
       -location_name => $hashref->{location_name},
       -adaptor => $self,
      );
  return $location;
}
