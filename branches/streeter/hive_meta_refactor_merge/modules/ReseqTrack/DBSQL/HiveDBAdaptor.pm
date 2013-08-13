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

  return  join(', ',qw(
      hive_db.hive_db_id
      hive_db.pipeline_id
      hive_db.name
      hive_db.host
      hive_db.port
      hive_db.created
      hive_db.retired
      hive_db.hive_version
      hive_db.is_seeded));
  
}

sub table_name{
  return "hive_db";
}


sub fetch_by_pipeline{

  my ($self, $pipeline, $forbid_retired) = @_;

  my $sql = "select ".$self->columns." from hive_db where pipeline_id = ?";
  if ($forbid_retired) {
    $sql .= " and retired is null";
  }

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
      . "(name, host, port, pipeline_id, hive_version, is_seeded,
          created) "
      . "values(?, ?, ?, ?, ?, ?, now()) ";
 
  my $sth = $self->prepare($sql);

  $sth->bind_param(1,  $hive_db->name);
  $sth->bind_param(2,  $hive_db->host);
  $sth->bind_param(3,  $hive_db->port);
  $sth->bind_param(4,  $hive_db->pipeline_id);
  $sth->bind_param(5,  $hive_db->hive_version);
  $sth->bind_param(6,  $hive_db->is_seeded // 0);
  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
 
  $sth->finish();
  $hive_db->dbID($dbID);
  $hive_db->adaptor($self);
  $hive_db->is_loaded(1);
}
#############


sub update{
  my ($self, $hive_db) = @_;

  my $sql = 
      "UPDATE hive_db SET name   = ?, ".
      "host     = ?, ".
      "port     = ?, ".
      "pipeline_id     = ?, ".
      "retired = ?, ".
      "is_seeded = ?, ".
      "hive_version = ?  ".
      "WHERE hive_db_id  = ?;";
  
  
  my $sth = $self->prepare($sql);
  
  
  $sth->bind_param(1,  $hive_db->name);
  $sth->bind_param(2,  $hive_db->host);
  $sth->bind_param(3,  $hive_db->port);
  $sth->bind_param(4,  $hive_db->pipeline_id);
  $sth->bind_param(5,  $hive_db->retired);
  $sth->bind_param(6,  $hive_db->is_seeded);
  $sth->bind_param(7,  $hive_db->hive_version);
  $sth->bind_param(8, $hive_db->dbID);
 
  $sth->execute();
  $sth->finish();
  $hive_db->is_loaded(1);
  return;
}

sub retire{
  my ($self, $hive_db) = @_;

  my $existing = $self->fetch_by_dbID($hive_db->dbID);
  throw('Does not exist in database: '. $hive_db->url)
      if !$existing;
  throw(join(' ', 'Already retired:', $hive_db->url, $hive_db->retired))
      if defined $hive_db->retired;

  my $sql = "select now()";
  my $sth = $self->prepare($sql);
 
  $sth->execute();
  $hive_db->retired($sth->fetchrow_arrayref->[0]);
  return $self->update($hive_db);
}

############


sub object_from_hashref{
    my ($self, $hashref) = @_;

    throw("Can't create an object from an undefined hashref") if(!$hashref);  

    my $hive_db = ReseqTrack::HiveDB->new
        (
         -dbID => $hashref->{hive_db_id},
         -adaptor => $self,
         -name             =>$hashref->{name},
         -port             =>$hashref->{port},
         -host             =>$hashref->{host},
         -pipeline_id      =>$hashref->{pipeline_id},
         -created          =>$hashref->{created},
         -retired          =>$hashref->{retired},
         -hive_version     =>$hashref->{hive_version},
         -is_seeded        =>$hashref->{is_seeded},
         -is_loaded            =>1,
     
    );
    return $hive_db; 
}


1;
