package ReseqTrack::DBSQL::HiveDBAdaptor;

use strict;
use warnings;

use base qw(ReseqTrack::DBSQL::BaseAdaptor);
use ReseqTrack::HiveDB;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument  qw(rearrange);


sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


sub columns{

  return  ' hive_db.hive_db_id, hive_db.url, hive_db.pipeline_id,  hive_db.created, hive_db.deleted, hive_db.hive_version ';
  
}

sub table_name{
  return "hive_db";
}



sub fetch_by_url{

  my ($self, $url) = @_;

  my $sql = "select ".$self->columns." from hive_db where url = ?";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $url);
  $sth->execute;
  my ($rowHashref) = $sth->fetchrow_hashref;
  my $hive_db = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;

  return $hive_db;
}

sub fetch_by_pipeline{

  my ($self, $pipeline) = @_;

  my $sql = "select ".$self->columns." from hive_db where pipeline_id = ?";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $pipeline->dbID);
  $sth->execute;
  my @hive_db;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $hive_db = $self->object_from_hashref($rowHashref);
    $hive_db->pipeline($pipeline);
    push(@hive_db, $hive_db);
  }
  $sth->finish;

  return \@hive_db;
}
#############


sub store{
  my ($self, $hive_db) = @_;
  my $sql = 
      "INSERT INTO  hive_db " 
      . "(url, pipeline_id, hive_version,
          created) "
      . "values(?, ?, ?, now()) ";
 
  my $sth = $self->prepare($sql);

  $sth->bind_param(1,  $hive_db->url);
  $sth->bind_param(2,  $hive_db->pipeline_id);
  $sth->bind_param(3,  $hive_db->hive_version);
  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
 
  $sth->finish();
  $hive_db->dbID($dbID);
  $hive_db->adaptor($self);
  $hive_db->loaded(1);
}
#############


sub update{
  my ($self, $hive_db) = @_;

  my $sql = 
      "UPDATE hive_db SET url   = ?, ".
      "pipeline_id     = ?, ".
      "deleted = ?, ".
      "hive_version = ?,  ".
      "WHERE pipeline_id  = ?;";
  
  
  my $sth = $self->prepare($sql);
  
  
  $sth->bind_param(1,  $hive_db->url);
  $sth->bind_param(2,  $hive_db->pipeline_id);
  $sth->bind_param(3,  $hive_db->deleted);
  $sth->bind_param(4,  $hive_db->hive_version);
  $sth->bind_param(5, $hive_db->dbID);
 
  $sth->execute();
  $sth->finish();
  $hive_db->loaded(1);
  return;
}
############


sub object_from_hashref{
    my ($self, $hashref) = @_;

    throw("Can't create an object from an undefined hashref") if(!$hashref);  

    my $hive_db = ReseqTrack::HiveDB->new
        (
         -dbID => $hashref->{hive_db_id},
         -adaptor => $self,
         -url              =>$hashref->{url},
         -pipeline_id      =>$hashref->{pipeline_id},
         -created          =>$hashref->{created},
         -deleted          =>$hashref->{deleted},
         -hive_version     =>$hashref->{hive_version},
         -loaded            =>1,
     
    );
    return $hive_db; 
}


1;
