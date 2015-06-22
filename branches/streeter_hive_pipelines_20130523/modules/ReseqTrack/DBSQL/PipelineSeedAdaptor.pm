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

  return  join(', ',qw(
      pipeline_seed.pipeline_seed_id
      pipeline_seed.seed_id
      pipeline_seed.hive_db_id
      pipeline_seed.is_running
      pipeline_seed.is_complete
      pipeline_seed.is_failed
      pipeline_seed.is_futile
      pipeline_seed.created
      pipeline_seed.completed));
  
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

  my $pipeline_seeds = $self->fetch_by_column_name('hive_db', $hive_db->dbID);
  foreach my $pipeline_seed (@$pipeline_seeds) {
    $pipeline_seed->hive_db($hive_db);
  }
  return $pipeline_seeds;
}

sub fetch_running_by_hive_db{
  my ($self, $hive_db) = @_;

  my $pipeline_seeds = $self->fetch_by_column_names(['hive_db_id', 'is_running'], [$hive_db->dbID, 1]);
  foreach my $pipeline_seed (@$pipeline_seeds) {
    $pipeline_seed->hive_db($hive_db);
  }
  return $pipeline_seeds;
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
  my $created_value = defined $pipeline_seed->created ? '?' : 'NOW()';
  my $sql = 
      "INSERT INTO  pipeline_seed " 
      . "(seed_id, hive_db_id, is_running, is_complete, is_failed, is_futile, completed, created) "
      . "values(?, ?, ?, ?, ?, ?, ?, $created_value) ";
 
  my $sth = $self->prepare($sql);

  $sth->bind_param(1,  $pipeline_seed->seed_id);
  $sth->bind_param(2,  $pipeline_seed->hive_db_id);
  $sth->bind_param(3,  $pipeline_seed->is_running);
  $sth->bind_param(4,  $pipeline_seed->is_complete // 0);
  $sth->bind_param(5,  $pipeline_seed->is_failed // 0);
  $sth->bind_param(6,  $pipeline_seed->is_futile // 0);
  $sth->bind_param(7,  $pipeline_seed->completed);
  $sth->bind_param(8,  $pipeline_seed->created) if defined $pipeline_seed->created;
  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
 
  $sth->finish();
  $pipeline_seed->dbID($dbID);
  $pipeline_seed->adaptor($self);
  $pipeline_seed->is_loaded(1);
}
#############


sub update{
  my ($self, $pipeline_seed) = @_;

  my $sql = 
      "UPDATE pipeline_seed SET ".
      "seed_id     = ?, ".
      "hive_db_id = ?, ".
      "is_running = ?,  ".
      "is_complete = ?,  ".
      "is_failed = ?,  ".
      "is_futile = ?,  ".
      "created = ?,  ".
      "completed = ?  ".
      "WHERE pipeline_seed_id  = ?;";
  
  
  my $sth = $self->prepare($sql);
  
  
  $sth->bind_param(1,  $pipeline_seed->seed_id);
  $sth->bind_param(2,  $pipeline_seed->hive_db_id);
  $sth->bind_param(3,  $pipeline_seed->is_running);
  $sth->bind_param(4,  $pipeline_seed->is_complete);
  $sth->bind_param(5,  $pipeline_seed->is_failed);
  $sth->bind_param(6,  $pipeline_seed->is_futile);
  $sth->bind_param(7,  $pipeline_seed->created);
  $sth->bind_param(8,  $pipeline_seed->completed);
  $sth->bind_param(9, $pipeline_seed->dbID);
  $pipeline_seed->is_loaded(1);
 
  $sth->execute();
  $sth->finish();
  return;
}

sub update_completed {
  my ($self, $pipeline_seed) = @_;
  $pipeline_seed->is_running(0);
  $pipeline_seed->is_complete(1);
  $pipeline_seed->is_failed(0);
  $pipeline_seed->is_futile(0);

  my $sql = 
      "UPDATE pipeline_seed SET ".
      "is_running = 0,  ".
      "is_complete = 1,  ".
      "is_failed = 0,  ".
      "is_futile = 0,  ".
      "completed = now()  ".
      "WHERE pipeline_seed_id  = ?;";
  
  my $sth = $self->prepare($sql);
  
  $sth->bind_param(1, $pipeline_seed->dbID);
 
  $sth->execute();
  $sth->finish();
  return;
}

sub update_failed {
  my ($self, $pipeline_seed, $is_futile) = @_;
  throw("must specify if failure is futile") if ! defined $is_futile;
  $pipeline_seed->is_running(0);
  $pipeline_seed->is_complete(0);
  $pipeline_seed->is_failed(1);
  $pipeline_seed->is_futile($is_futile);

  my $sql = 
      "UPDATE pipeline_seed SET ".
      "is_running = 0,  ".
      "is_complete = 0,  ".
      "is_failed = 1,  ".
      "is_futile = ?,  ".
      "completed = now()  ".
      "WHERE pipeline_seed_id  = ?;";
  
  my $sth = $self->prepare($sql);
  
  $sth->bind_param(1, $pipeline_seed->dbID);
  $sth->bind_param(2, $is_futile);
 
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
         -is_running          =>$hashref->{is_running},
         -is_complete          =>$hashref->{is_complete},
         -is_failed          =>$hashref->{is_failed},
         -is_futile          =>$hashref->{is_futile},
         -created          =>$hashref->{created},
         -completed          =>$hashref->{completed},
         -is_loaded            =>1,
     
    );
    return $pipeline_seed; 
}


1;
