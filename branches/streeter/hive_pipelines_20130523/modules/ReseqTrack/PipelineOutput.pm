package ReseqTrack::PipelineOutput;

use strict;
use warnings;

use base qw(ReseqTrack::Base);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);



sub new{
  my ($class, @args) = @_;
  my  $self = $class->SUPER::new(@args);
  my ($pipeline_seed, $pipeline_seed_id,
      $hive_db, $hive_db_id,
      $pipeline, $pipeline_id,
      $seed, $seed_id,
      $output, $output_id,
      $action, $table_name, $loaded,
      ) = 
      rearrange(['PIPELINE_SEED', 'PIPELINE_SEED_ID',
                'HIVE_DB', 'HIVE_DB_ID',
                'PIPELINE', 'PIPELINE_ID',
                'SEED', 'SEED_ID',
                'OUTPUT', 'OUTPUT_ID',
                'ACTION', 'TABLE_NAME', 'LOADED',
                 ], @args);
  
  #don't trigger any lazy loading if parameters are undefined
  $self->pipeline_seed($pipeline_seed) if $pipeline_seed;
  $self->pipeline_seed_id($pipeline_seed_id) if $pipeline_seed_id;
  $self->hive_db($hive_db) if $hive_db;
  $self->hive_db_id($hive_db_id) if $hive_db_id;
  $self->pipeline($pipeline) if $pipeline;
  $self->pipeline_id($pipeline_id) if $pipeline_id;
  $self->seed($seed) if $seed;
  $self->seed_id($seed_id) if $seed_id;
  $self->output($output) if $output;
  $self->output_id($output_id) if $output;
  $self->action($action) if $action;
  $self->table_name($table_name) if $table_name;
  $self->loaded($loaded);

  return $self;
}


# Controls lazy loading
sub loaded{
  my ($self, $loaded) = @_;
  if (@_>1){
    $self->{loaded} = $loaded || 0;
  }
  return $self->{loaded};
}

sub load {
  my ($self) = @_;
  return if !$self->adaptor;
  my $object = $self->adaptor->fetch_by_dbID($self->dbID);
  $self->loaded(1); # set to 1 now to prevent recursive loops
  $self->table_name($self->table_name // $object->table_name);
  $self->pipeline_seed_id($self->pipeline_seed_id // $object->pipeline_seed_id);
  $self->output_id($self->output_id // $object->output_id);
  $self->action($self->action // $object->action);

}

sub pipeline_seed{
  my ($self, $pipeline_seed) = @_;
  if($pipeline_seed){
    throw("not a ReseqTrack::PipelineSeed object") if ref($pipeline_seed) ne 'ReseqTrack::PipelineSeed';
    $self->{pipeline_seed} = $pipeline_seed;
  }
  if (!$self->{pipeline_seed} && !$self->loaded) {
    $self->load;
  }
  if (!$self->{pipeline_seed}) {
    my $adaptor = $self->adaptor ? $self->adaptor->db->get_PipelineSeedAdaptor : undef;
    $self->{pipeline_seed} = ReseqTrack::PipelineSeed->new(-adaptor => $adaptor);
  }
  return $self->{pipeline_seed};
}

sub pipeline_seed_id{
  my ($self, $pipeline_seed_id) = @_;
  if($pipeline_seed_id){
    if ($self->{pipeline_seed}) {
      $self->pipeline_seed->dbID($pipeline_seed_id);
    }
    else {
      my $adaptor = $self->adaptor ? $self->adaptor->db->get_PipelineSeedAdaptor : undef;
      my $pipeline_seed = ReseqTrack::PipelineSeed->new(-dbID => $pipeline_seed_id, -adaptor => $adaptor);
      $self->pipeline_seed($pipeline_seed);
    }
  }
  return $self->pipeline_seed->dbID;
}

sub hive_db{
  my ($self, $hive_db) = @_;
  return $self->pipeline_seed->hive_db($hive_db);
}

sub hive_db_id{
  my ($self, $hive_db_id) = @_;
  return $self->pipeline_seed->hive_db_id($hive_db_id);
}

sub pipeline{
  my ($self, $pipeline) = @_;
  return $self->hive_db->pipeline($pipeline);
}

sub pipeline_id{
  my ($self, $pipeline_id) = @_;
  return $self->hive_db->pipeline_id($pipeline_id);
}

sub seed{
  my ($self, $seed) = @_;
  return $self->pipeline_seed->seed($seed);
}

sub seed_id{
  my ($self, $seed_id) = @_;
  return $self->pipeline_seed->seed_id($seed_id);
}



sub output{
  my ($self, $output) = @_;
  if($output){
    throw("not a ReseqTrack object") if ref($output) !~ /^ReseqTrack::/;
    $self->{output} = $output;
  }
  if (!$self->{output}) {
    my $table_name = $self->table_name;
    my $adaptor = $self->adaptor;
    my $output_id = $self->output_id;
    if ($table_name && $adaptor && $output_id) {
      my $output = $table_name eq 'file' ? $adaptor->db->get_FileAdaptor->fetch_by_dbID($output_id)
               : $table_name eq 'collection' ? $adaptor->db->get_CollectionAdaptor->fetch_by_dbID($output_id)
               : $table_name eq 'run_meta_info' ? $adaptor->db->get_RunMetaInfoAdaptor->fetch_by_dbID($output_id)
               : $table_name eq 'input_string' ? $adaptor->db->get_InputStringAdaptor->fetch_by_dbID($output_id)
               : throw("cannot get seed for table $table_name");
      $self->{output} = $output;
    }
  }
  return $self->{output};
}

sub output_id{
  my ($self, $output_id) = @_;
  if($output_id){
    $self->{output_id} = $output_id;
    if ($self->{output} && $self->output->dbID != $output_id) {
      $self->{output} = undef;
    }
  }
  if (!$self->{output_id} && $self->{output}) {
    $self->{output_id} = $self->output->dbID;
  }
  if (!$self->{output_id} && !$self->loaded) {
    $self->load($self);
  }
  return $self->{output_id};
}


sub table_name{
  my ($self, $table_name) = @_;
  if($table_name){
    $self->{table_name} = $table_name;
  }
  if (!$self->{table_name} && $self->{output}) {
    $self->{table_name} = $self->output->object_table_name;
  }
  if (!$self->{table_name} && !$self->loaded) {
    $self->load($self);
  }
  return $self->{table_name};
}

sub action{
  my ($self, $action) = @_;
  if($action){
    $self->{action} = $action;
  }
  if (!$self->{action} && !$self->loaded) {
    $self->load($self);
  }
  return $self->{action};
}


sub object_table_name{
  my ($self) = @_;
  return 'pipeline_output';
}

1;

