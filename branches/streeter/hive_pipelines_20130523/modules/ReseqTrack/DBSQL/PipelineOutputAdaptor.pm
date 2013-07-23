package ReseqTrack::DBSQL::PipelineOutputAdaptor;

use strict;
use warnings;

use base qw(ReseqTrack::DBSQL::BaseAdaptor);
use ReseqTrack::PipelineOutput;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument  qw(rearrange);


sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


sub columns{

  return  ' pipeline_output.pipeline_output_id, pipeline_output.pipeline_seed_id, pipeline_output.table_name, pipeline_output.output_id, pipeline_output.action ';
  
}

sub table_name{
  return "pipeline_output";
}

sub fetch_by_pipeline_seed{

  my ($self, $pipeline_seed) = @_;

  my $sql = "select ".$self->columns." from pipeline_output where pipeline_seed_id = ? ";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $pipeline_seed->dbID);
  $sth->execute;
  my @pipeline_outputs;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $pipeline_output = $self->object_from_hashref($rowHashref);
    $pipeline_output->pipeline_seed($pipeline_seed);
    push(@pipeline_outputs, $pipeline_output);
  }
  $sth->finish;

  return \@pipeline_outputs;
}

sub fetch_by_output{

  my ($self, $output) = @_;

  my $sql = "select ".$self->columns." from pipeline_output where table_name = ? and output_id = ? ";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $output->object_table_name);
  $sth->bind_param(2, $output->dbID);
  $sth->execute;
  my @pipeline_outputs;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $pipeline_output = $self->object_from_hashref($rowHashref);
    $pipeline_output->output($output);
    push(@pipeline_outputs, $pipeline_output);
  }
  $sth->finish;

  return \@pipeline_outputs;
}

sub fetch_by_hive_db{

  my ($self, $hive_db) = @_;

  my $sql = "select ".$self->columns." from pipeline_output, pipeline_seed"
    ." where pipeline_output.pipeline_seed_id = pipeline_seed.pipeline_seed_id"
    ." and pipeline_seed.hive_db_id= ? ";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $hive_db->dbID);
  $sth->execute;
  my @pipeline_outputs;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $pipeline_output = $self->object_from_hashref($rowHashref);
    $pipeline_output->hive_db($hive_db);
    push(@pipeline_outputs, $pipeline_output);
  }
  $sth->finish;

  return \@pipeline_outputs;
}

sub fetch_by_pipeline{

  my ($self, $pipeline) = @_;

  my $sql = "select ".$self->columns.", hive_db.hive_db_id from pipeline_output, pipeline_seed, hive_db"
    ." where pipeline_output.pipeline_seed_id = pipeline_seed.pipeline_seed_id"
    ." and pipeline_seed.hive_db_id = hive_db.hive_db_id"
    ." and hive_db.pipeline_id= ? ";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $pipeline->dbID);
  $sth->execute;
  my @pipeline_outputs;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $pipeline_output = $self->object_from_hashref($rowHashref);
    $pipeline_output->hive_db_id($rowHashref->{hive_db_id});
    $pipeline_output->pipeline($pipeline);
    push(@pipeline_outputs, $pipeline_output);
  }
  $sth->finish;

  return \@pipeline_outputs;
}




#############


sub store{
  my ($self, $pipeline_output) = @_;
  my $sql = 
      "INSERT INTO  pipeline_output " 
      . "(pipeline_seed_id, table_name, output_id, action) "
      . "values(?, ?, ?, ?) ";
 
  my $sth = $self->prepare($sql);

  $sth->bind_param(1,  $pipeline_output->pipeline_seed_id);
  $sth->bind_param(2,  $pipeline_output->table_name);
  $sth->bind_param(3,  $pipeline_output->output_id);
  $sth->bind_param(4,  $pipeline_output->action);
  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
 
  $sth->finish();
  $pipeline_output->dbID($dbID);
  $pipeline_output->adaptor($self);
  $pipeline_output->loaded(1);
}
#############


sub update{
  my ($self, $pipeline_output) = @_;

  my $sql = 
      "UPDATE pipeline_output SET ".
      "pipeline_seed_id     = ?, ".
      "table_name = ?, ".
      "output_id = ?,  ".
      "action = ?,  ".
      "WHERE pipeline_output_id  = ?;";
  
  
  my $sth = $self->prepare($sql);
  
  
  $sth->bind_param(1,  $pipeline_output->pipeline_seed_id);
  $sth->bind_param(2,  $pipeline_output->table_name);
  $sth->bind_param(3,  $pipeline_output->output_id);
  $sth->bind_param(4,  $pipeline_output->action);
  $sth->bind_param(5, $pipeline_output->dbID);
 
  $sth->execute();
  $sth->finish();
  $pipeline_output->loaded(1);
  return;
}
############


sub object_from_hashref{
    my ($self, $hashref) = @_;

    throw("Can't create an object from an undefined hashref") if(!$hashref);  

    my $pipeline_output = ReseqTrack::PipelineOutput->new
        (
         -dbID => $hashref->{pipeline_output_id},
         -adaptor => $self,
         -pipeline_seed_id =>$hashref->{pipeline_seed_id},
         -output_id       =>$hashref->{output_id},
         -table_name          =>$hashref->{table_name},
         -action          =>$hashref->{action},
         -loaded            =>1,
     
    );
    return $pipeline_output; 
}


1;
