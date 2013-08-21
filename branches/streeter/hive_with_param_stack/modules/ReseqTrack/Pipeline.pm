package ReseqTrack::Pipeline;

use strict;
use warnings;

use base qw(ReseqTrack::HasHistory);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);



sub new{
  my ($class, @args) = @_;
  my  $self = $class->SUPER::new(@args);
  my ($name, $table_name,
      $config_module, $config_options,
      $seeding_module, $seeding_options,
      $created, $is_loaded) = 
      rearrange(['NAME', 'TABLE_NAME',
                 'CONFIG_MODULE', 'CONFIG_OPTIONS',
                 'SEEDING_MODULE', 'SEEDING_OPTIONS',
                 'CREATED', 'IS_LOADED'], @args);
  
  $self->name($name) if $name;
  $self->table_name($table_name) if $table_name;
  $self->config_module($config_module) if $config_module;
  $self->config_options($config_options) if $config_options;
  $self->seeding_module($seeding_module) if $seeding_module;
  $self->seeding_options($seeding_options) if $seeding_options;
  $self->created($created) if $created;
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
  $self->is_loaded(1);
  $self->name($self->name // $object->name);
  $self->table_name($self->table_name // $object->table_name);
  $self->config_module($self->config_module // $object->config_module);
  $self->config_options($self->config_options // $object->config_options);
  $self->seeding_module($self->seeding_module // $object->seeding_module);
  $self->seeding_options($self->seeding_options // $object->seeding_options);
  $self->created($self->created // $object->created);
}

sub name{
  my ($self, $name) = @_;
  if($name){
    $self->{name} = $name;
  }
  if (!$self->{name} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{name};
}

sub table_name{
  my ($self, $table_name) = @_;
  if($table_name){
    $self->{table_name} = $table_name;
  }
  if (!$self->{table_name} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{table_name};
}

sub config_module{
  my ($self, $config_module) = @_;
  if($config_module){
    $self->{config_module} = $config_module;
  }
  if (!$self->{config_module} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{config_module};
}

sub config_options{
  my ($self, $config_options) = @_;
  if($config_options){
    $self->{config_options} = $config_options;
  }
  if (!$self->{config_options} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{config_options};
}

sub seeding_module{
  my ($self, $seeding_module) = @_;
  if($seeding_module){
    $self->{seeding_module} = $seeding_module;
  }
  if (!$self->{seeding_module} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{seeding_module};
}

sub seeding_options{
  my ($self, $seeding_options) = @_;
  if($seeding_options){
    $self->{seeding_options} = $seeding_options;
  }
  if (!$self->{seeding_options} && !$self->is_loaded) {
    $self->load;
  }
  return $self->{seeding_options};
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


sub object_table_name{
  my ($self) = @_;
  return 'pipeline';
}

1;

