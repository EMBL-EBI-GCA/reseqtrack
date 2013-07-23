package ReseqTrack::DBSQL::PipelineSeedAdaptor;

use strict;
use warnings;

use base qw(ReseqTrack::DBSQL::BaseAdaptor);
use ReseqTrack::PipelineSeed;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument  qw(rearrange);


sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


sub columns{

  return  ' pipeline_seed.pipeline_seed_id, pipeline_seed.seed_id, pipeline_seed.hive_db_id, pipeline_seed.status, pipeline_seed.updated ';
  
}

sub table_name{
  return "pipeline_seed";
}

sub fetch_by_seed{

  my ($self, $seed) = @_;

  my $sql = "select ".$self->columns.", hive_db.pipeline_id from pipeline_seed, hive_db, pipeline"
          ." where pipeline_seed.hive_db_id = hive_db.hive_db_id"
          ." and hive_db.pipeline_id = pipeline.pipeline_id"
          ." and pipeline.table_name = ? and seed_id = ?";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $seed->object_table_name);
  $sth->bind_param(2, $seed->dbID);
  $sth->execute;
  my @pipeline_seeds;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $pipeline_seed = $self->object_from_hashref($rowHashref);
    $pipeline_seed->seed($seed);
    $pipeline_seed->pipeline_id($rowHashref->{pipeline_id});
    push(@pipeline_seeds, $pipeline_seed);
  }
  $sth->finish;

  return \@pipeline_seeds;
}

sub fetch_by_hive_db{

  my ($self, $hive_db) = @_;

  my $sql = "select ".$self->columns." from pipeline_seed where hive_db_id = ?";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $hive_db->dbID);
  $sth->execute;
  my @pipeline_seeds;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $pipeline_seed = $self->object_from_hashref($rowHashref);
    $pipeline_seed->hive_db($hive_db);
    push(@pipeline_seeds, $pipeline_seed);
  }
  $sth->finish;

  return \@pipeline_seeds;
}

sub fetch_by_pipeline{

  my ($self, $pipeline) = @_;

  my $sql = "select ".$self->columns." from pipeline_seed, hive_db where pipeline_seed.hive_db_id = hive_db.hive_db_id and hive_db.pipeline_id = ?";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $pipeline->dbID);
  $sth->execute;
  my @pipeline_seeds;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $pipeline_seed = $self->object_from_hashref($rowHashref);
    $pipeline_seed->pipeline($pipeline);
    push(@pipeline_seeds, $pipeline_seed);
  }
  $sth->finish;

  return \@pipeline_seeds;
}

#############


sub store{
  my ($self, $pipeline_seed) = @_;
  my $sql = 
      "INSERT INTO  pipeline_seed " 
      . "(seed_id, hive_db_id, status, updated) "
      . "values(?, ?, ?, now()) ";
 
  my $sth = $self->prepare($sql);

  $sth->bind_param(1,  $pipeline_seed->seed_id);
  $sth->bind_param(2,  $pipeline_seed->hive_db_id);
  $sth->bind_param(3,  $pipeline_seed->status);
  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
 
  $sth->finish();
  $pipeline_seed->dbID($dbID);
  $pipeline_seed->adaptor($self);
  $pipeline_seed->loaded(1);
}
#############


sub update{
  my ($self, $pipeline_seed) = @_;

  my $sql = 
      "UPDATE pipeline_seed SET ".
      "seed_id     = ?, ".
      "hive_db_id = ?, ".
      "status = ?,  ".
      "updated = now(),  ".
      "WHERE pipeline_seed_id  = ?;";
  
  
  my $sth = $self->prepare($sql);
  
  
  $sth->bind_param(1,  $pipeline_seed->seed_id);
  $sth->bind_param(2,  $pipeline_seed->hive_db_id);
  $sth->bind_param(3,  $pipeline_seed->status);
  $sth->bind_param(4, $pipeline_seed->dbID);
  $pipeline_seed->loaded(1);
 
  $sth->execute();
  $sth->finish();
  return;
}
############


sub object_from_hashref{
    my ($self, $hashref) = @_;

    throw("Can't create an object from an undefined hashref") if(!$hashref);  

    my $pipeline_seed = ReseqTrack::PipelineSeed->new
        (
         -dbID => $hashref->{pipeline_seed_id},
         -adaptor => $self,
         -seed_id              =>$hashref->{seed_id},
         -hive_db_id         =>$hashref->{hive_db_id},
         -status          =>$hashref->{status},
         -updated          =>$hashref->{updated},
         -loaded            =>1,
     
    );
    return $pipeline_seed; 
}


1;
