package ReseqTrack::PipelineSeed;

use strict;
use warnings;

use base qw(ReseqTrack::HasHistory);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);



sub new{
  my ($class, @args) = @_;
  my  $self = $class->SUPER::new(@args);
  my ($seed, $seed_id,
      $hive_db, $hive_db_id,
      $pipeline, $pipeline_id,
      $is_running, $is_complete, $is_failed, $is_futile,
      $created, $completed, $is_loaded) = 
      rearrange(['SEED', 'SEED_ID',
                'HIVE_DB', 'HIVE_DB_ID',
                'PIPELINE', 'PIPELINE_ID',
                'IS_RUNNING', 'IS_COMPLETE', 'IS_FAILED', 'IS_FUTILE',
                'CREATED', 'COMPLETED', 'IS_LOADED',
                 ], @args);
  
  #don't trigger any lazy loading if parameters are undefined
  $self->seed($seed) if $seed;
  $self->seed_id($seed_id) if $seed_id;
  $self->hive_db($hive_db) if $hive_db;
  $self->hive_db_id($hive_db_id) if $hive_db_id;
  $self->pipeline($pipeline) if $pipeline;
  $self->pipeline_id($pipeline_id) if $pipeline_id;
  $self->is_running($is_running) if defined $is_running;
  $self->is_complete($is_complete) if defined $is_complete;
  $self->is_failed($is_failed) if defined $is_failed;
  $self->is_futile($is_futile) if defined $is_futile;
  $self->created($created) if $created;
  $self->completed($completed) if $completed;
  $self->is_loaded($is_loaded);

  return $self;
}

# Controls lazy loading
sub is_loaded{
  my ($self, $is_loaded) = @_;
  if (@_>1){
    $self->{is_loaded} = $is_loaded || 0;
  }
  return $self->{is_loaded};
}

sub load {
  my ($self) = @_;
  return if !$self->adaptor;
  my $object = $self->adaptor->fetch_by_dbID($self->dbID);
  $self->is_loaded(1); # set to 1 now to prevent recursive loops
  $self->seed_id($self->seed_id // $object->seed_id);
  $self->hive_db_id($self->hive_db_id // $object->hive_db_id);
  $self->is_running($self->is_running // $object->is_running);
  $self->is_complete($self->is_complete // $object->is_complete);
  $self->is_failed($self->is_failed // $object->is_failed);
  $self->is_futile($self->is_futile // $object->is_futile);
  $self->created($self->created // $object->created);
  $self->completed($self->completed // $object->completed);

}


sub hive_db{
  my ($self, $hive_db) = @_;
  if($hive_db){
    throw("not a ReseqTrack::HiveDB object") if ref($hive_db) ne 'ReseqTrack::HiveDB';
    $self->{hive_db} = $hive_db;
  }
  if (!$self->{hive_db} && !$self->is_loaded) {
    $self->load;
  }
  if (!$self->{hive_db}) {
    my $adaptor = $self->adaptor ? $self->adaptor->db->get_HiveDBAdaptor : undef;
    $self->{hive_db} = ReseqTrack::HiveDB->new(-adaptor => $adaptor);
  }
  return $self->{hive_db};
}

sub hive_db_id{
  my ($self, $hive_db_id) = @_;
  if($hive_db_id){
    if ($self->{hive_db}) {
      $self->hive_db->dbID($hive_db_id);
    }
    else {
      my $adaptor = $self->adaptor ? $self->adaptor->db->get_HiveDBAdaptor : undef;
      my $hive_db = ReseqTrack::HiveDB->new(-dbID => $hive_db_id, -adaptor => $adaptor);
      $self->hive_db($hive_db);
    }
  }
  return $self->hive_db->dbID;
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
  if($seed){
    throw("not a ReseqTrack object") if ref($seed) !~ /^ReseqTrack::/;
    $self->{seed} = $seed;
  }
  if (!$self->{seed}) {
    my $table_name = $self->pipeline->table_name;
    my $adaptor = $self->adaptor;
    my $seed_id = $self->seed_id;
    if ($table_name && $adaptor && $seed_id) {
      my $seed = $table_name eq 'file' ? $adaptor->db->get_FileAdaptor->fetch_by_dbID($seed_id)
               : $table_name eq 'collection' ? $adaptor->db->get_CollectionAdaptor->fetch_by_dbID($seed_id)
               : $table_name eq 'run_meta_info' ? $adaptor->db->get_RunMetaInfoAdaptor->fetch_by_dbID($seed_id)
               : $table_name eq 'input_string' ? $adaptor->db->get_InputStringAdaptor->fetch_by_dbID($seed_id)
               : throw("cannot get seed for table $table_name");
      $self->{seed} = $seed;
    }
  }
  return $self->{seed};
}

sub seed_id{
  my ($self, $seed_id) = @_;
  if($seed_id){
    $self->{seed_id} = $seed_id;
    if ($self->{seed} && $self->seed->dbID != $seed_id) {
      $self->{seed} = undef;
    }
  }
  if (!$self->{seed_id} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{seed_id};
}

sub is_running{
  my ($self, $is_running) = @_;
  if(defined $is_running){
    $self->{is_running} = $is_running;
  }
  if (! defined $self->{is_running} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{is_running};
}

sub is_complete{
  my ($self, $is_complete) = @_;
  if(defined $is_complete){
    $self->{is_complete} = $is_complete;
  }
  if (! defined $self->{is_complete} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{is_complete};
}

sub is_failed{
  my ($self, $is_failed) = @_;
  if(defined $is_failed){
    $self->{is_failed} = $is_failed;
  }
  if (! defined $self->{is_failed} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{is_failed};
}

sub is_futile{
  my ($self, $is_futile) = @_;
  if(defined $is_futile){
    $self->{is_futile} = $is_futile;
  }
  if (! defined $self->{is_futile} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{is_futile};
}

sub created{
  my ($self, $created) = @_;
  if($created){
    $self->{created} = $created;
  }
  if (!$self->{created} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{created};
}

sub completed{
  my ($self, $completed) = @_;
  if($completed){
    $self->{completed} = $completed;
  }
  if (!$self->{completed} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{completed};
}


sub object_table_name{
  my ($self) = @_;
  return 'pipeline_seed';
}

1;

