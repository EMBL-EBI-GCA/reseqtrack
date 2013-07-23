package ReseqTrack::PipelineSeed;

use strict;
use warnings;

use base qw(ReseqTrack::Base);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);



sub new{
  my ($class, @args) = @_;
  my  $self = $class->SUPER::new(@args);
  my ($seed, $seed_id,
      $hive_db, $hive_db_id,
      $pipeline, $pipeline_id,
      $status, $updated, $loaded) = 
      rearrange(['SEED', 'SEED_ID',
                'HIVE_DB', 'HIVE_DB_ID',
                'PIPELINE', 'PIPELINE_ID',
                'STATUS', 'UPDATED', 'LOADED',
                 ], @args);
  
  #don't trigger any lazy loading if parameters are undefined
  $self->seed($seed) if $seed;
  $self->seed_id($seed_id) if $seed_id;
  $self->hive_db($hive_db) if $hive_db;
  $self->hive_db_id($hive_db_id) if $hive_db_id;
  $self->pipeline($pipeline) if $pipeline;
  $self->pipeline_id($pipeline_id) if $pipeline_id;
  $self->status($status) if $status;
  $self->updated($updated) if $updated;
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
  $self->seed_id($self->seed_id // $object->seed_id);
  $self->hive_db_id($self->hive_db_id // $object->hive_db_id);
  $self->status($self->status // $object->status);
  $self->updated($self->updated // $object->updated);

}


sub hive_db{
  my ($self, $hive_db) = @_;
  if($hive_db){
    throw("not a ReseqTrack::HiveDB object") if ref($hive_db) ne 'ReseqTrack::HiveDB';
    $self->{hive_db} = $hive_db;
  }
  if (!$self->{hive_db} && !$self->loaded) {
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
    print $table_name, "\n";
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
  if (!$self->{seed_id} && !$self->loaded) {
    $self->adaptor->load($self);
  }
  return $self->{seed_id};
}

sub status{
  my ($self, $status) = @_;
  if($status){
    $self->{status} = $status;
  }
  if (!$self->{status} && !$self->loaded) {
    $self->adaptor->load($self);
  }
  return $self->{status};
}

sub updated{
  my ($self, $updated) = @_;
  if($updated){
    $self->{updated} = $updated;
  }
  if (!$self->{updated} && !$self->loaded) {
    $self->adaptor->load($self);
  }
  return $self->{updated};
}


sub object_table_name{
  my ($self) = @_;
  return 'pipeline_seed';
}

1;

