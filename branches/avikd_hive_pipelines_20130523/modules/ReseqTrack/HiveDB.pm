package ReseqTrack::HiveDB;

use strict;
use warnings;

use base qw(ReseqTrack::Base);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);


sub new{
  my ($class, @args) = @_;
  my  $self = $class->SUPER::new(@args);
  my ($url, $pipeline_id, $pipeline, $created, $retired, $hive_version, $is_seeded, $is_loaded,
      ) = 
      rearrange(['URL', 'PIPELINE_ID', 'PIPELINE', 'CREATED', 'RETIRED', 'HIVE_VERSION', 'IS_SEEDED', 'IS_LOADED',
                 ], @args);
  
  #don't trigger any lazy loading if parameters are undefined
  $self->url($url) if $url;
  $self->pipeline($pipeline) if $pipeline;
  $self->pipeline_id($pipeline_id) if $pipeline_id;
  $self->created($created) if $created;
  $self->retired($retired) if $retired;
  $self->hive_version($hive_version) if $hive_version;
  $self->is_seeded($is_seeded) if defined $is_seeded;
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
  $self->url($self->url // $object->url);
  $self->pipeline_id($self->pipeline_id // $object->pipeline_id);
  $self->created($self->created // $object->created);
  $self->retired($self->retired // $object->retired);
  $self->hive_version($self->hive_version // $object->hive_version);
  $self->is_seeded($self->is_seeded // $object->is_seeded);
}


sub url{
  my ($self, $url) = @_;
  if($url){
    $self->{url} = $url;
  }
  if (!$self->{url} && !$self->is_loaded) {
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
  if (!$self->{pipeline} && !$self->is_loaded) {
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
  if (!$self->{created} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{created};
}

sub retired{
  my ($self, $retired) = @_;
  if($retired){
    $self->{retired} = $retired;
  }
  if (!$self->{retired} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{retired};
}

sub hive_version{
  my ($self, $hive_version) = @_;
  if($hive_version){
    $self->{hive_version} = $hive_version;
  }
  if (!$self->{hive_version} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{hive_version};
}

sub is_seeded{
  my ($self, $is_seeded) = @_;
  if(defined $is_seeded){
    $self->{is_seeded} = $is_seeded;
  }
  if (!defined $self->{is_seeded} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{is_seeded};
}



sub object_table_name{
  my ($self) = @_;
  return 'hive_db';
}

1;

