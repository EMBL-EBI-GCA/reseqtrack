package ReseqTrack::HiveDB;

use strict;
use warnings;

use base qw(ReseqTrack::Base);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);



sub new{
  my ($class, @args) = @_;
  my  $self = $class->SUPER::new(@args);
  my ($url, $pipeline_id, $pipeline, $created, $deleted, $hive_version, $loaded,
      ) = 
      rearrange(['URL', 'PIPELINE_ID', 'PIPELINE', 'CREATED', 'DELETED', 'HIVE_VERSION', 'LOADED',
                 ], @args);
  
  #don't trigger any lazy loading if parameters are undefined
  $self->url($url) if $url;
  $self->pipeline($pipeline) if $pipeline;
  $self->pipeline_id($pipeline_id) if $pipeline_id;
  $self->created($created) if $created;
  $self->deleted($deleted) if $deleted;
  $self->hive_version($hive_version) if $hive_version;
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
  $self->url($self->url // $object->url);
  $self->pipeline_id($self->pipeline_id // $object->pipeline_id);
  $self->created($self->created // $object->created);
  $self->deleted($self->deleted // $object->deleted);
  $self->hive_version($self->hive_version // $object->hive_version);
}


sub url{
  my ($self, $url) = @_;
  if($url){
    $self->{url} = $url;
  }
  if (!$self->{url} && !$self->loaded) {
    $self->load;
  }
  return $self->{url};
}

sub pipeline{
  my ($self, $pipeline) = @_;
  if($pipeline){
    throw("not a ReseqTrack::Pipeline object") if ref($pipeline) ne 'ReseqTrack::Pipeline';
    $self->{pipeline} = $pipeline;
  }
  if (!$self->{pipeline} && !$self->loaded) {
    $self->load;
  }
  if (!$self->{pipeline}) {
    my $adaptor = $self->adaptor ? $self->adaptor->db->get_PipelineAdaptor : undef;
    $self->{pipeline} = ReseqTrack::Pipeline->new(-adaptor => $adaptor);
  }
  return $self->{pipeline};
}

sub pipeline_id{
  my ($self, $pipeline_id) = @_;
  if($pipeline_id){
    if ($self->{pipeline}) {
      $self->pipeline->dbID($pipeline_id);
    }
    else {
      my $adaptor = $self->adaptor ? $self->adaptor->db->get_PipelineAdaptor : undef;
      my $pipeline = ReseqTrack::Pipeline->new(-dbID => $pipeline_id, -adaptor => $adaptor);
      $self->pipeline($pipeline);
    }
  }
  return $self->pipeline->dbID;
}


sub created{
  my ($self, $created) = @_;
  if($created){
    $self->{created} = $created;
  }
  if (!$self->{created} && !$self->loaded) {
    $self->adaptor->load($self);
  }
  return $self->{created};
}

sub deleted{
  my ($self, $deleted) = @_;
  if($deleted){
    $self->{deleted} = $deleted;
  }
  if (!$self->{deleted} && !$self->loaded) {
    $self->adaptor->load($self);
  }
  return $self->{deleted};
}

sub hive_version{
  my ($self, $hive_version) = @_;
  if($hive_version){
    $self->{hive_version} = $hive_version;
  }
  if (!$self->{hive_version} && !$self->loaded) {
    $self->adaptor->load($self);
  }
  return $self->{hive_version};
}



sub object_table_name{
  my ($self) = @_;
  return 'hive_db';
}

1;

