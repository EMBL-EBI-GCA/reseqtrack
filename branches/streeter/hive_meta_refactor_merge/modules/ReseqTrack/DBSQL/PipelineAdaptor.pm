package ReseqTrack::DBSQL::PipelineAdaptor;

use strict;
use warnings;

use base qw(ReseqTrack::DBSQL::BaseAdaptor);
use ReseqTrack::Pipeline;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument  qw(rearrange);


sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


sub columns{

  return  join(', ',qw(
      pipeline.pipeline_id
      pipeline.name
      pipeline.table_name
      pipeline.config_module
      pipeline.config_options
      pipeline.created
      ));
}

sub table_name{
  return "pipeline";
}


sub fetch_by_name{

  my ($self, $name) = @_;

  my $sql = "select ".$self->columns." from pipeline where name = ?";

  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name);
  $sth->execute;
  my ($rowHashref) = $sth->fetchrow_hashref;
  my $pipeline = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;

  return $pipeline;
}
#############


sub store{
  my ($self, $pipeline) = @_;
  my $sql = 
      'INSERT INTO  pipeline' 
      . ' (name, table_name,'
      . ' config_module, config_options,'
      . ' created)'
      . ' values(?, ?, ?, ?, ?, ?, now()) ';
 
  my $sth = $self->prepare($sql);

  $sth->bind_param(1,  $pipeline->name);
  $sth->bind_param(2,  $pipeline->table_name);
  $sth->bind_param(3,  $pipeline->config_module);
  $sth->bind_param(4,  $pipeline->config_options);
  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
 
  $sth->finish();
  $pipeline->dbID($dbID);
  $pipeline->adaptor($self);
  $self->store_history($pipeline);
  $pipeline->is_loaded(1);
}
#############


sub update{
  my ($self, $pipeline) = @_;

  my $sql = 
      "UPDATE pipeline SET name   = ?, ".
      "table_name     = ?, ".
      "config_module = ?,  ".
      "config_options = ?,  ".
      "WHERE pipeline_id  = ?";
  
  
  my $sth = $self->prepare($sql);
  
  
  $sth->bind_param(1,  $pipeline->name);
  $sth->bind_param(2,  $pipeline->table_name);
  $sth->bind_param(3,  $pipeline->config_module);
  $sth->bind_param(4,  $pipeline->config_options);
  $sth->bind_param(5, $pipeline->dbID);
 
  $sth->execute();
  $sth->finish();
  $self->store_history($pipeline);
  $pipeline->is_loaded(1);
  return;
}
############


sub object_from_hashref{
    my ($self, $hashref) = @_;

    throw("Can't create an object from an undefined hashref") if(!$hashref);  

    my $pipeline = ReseqTrack::Pipeline->new
        (
         -dbID => $hashref->{pipeline_id},
         -adaptor => $self,
         -name              =>$hashref->{name},
         -table_name        =>$hashref->{table_name},
         -config_module     =>$hashref->{config_module},
         -config_options    =>$hashref->{config_options},
         -created           =>$hashref->{created},
         -is_loaded         =>1,
     
    );
    return $pipeline; 
}


1;
