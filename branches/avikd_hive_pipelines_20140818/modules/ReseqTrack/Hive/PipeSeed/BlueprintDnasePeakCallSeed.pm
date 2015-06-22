package ReseqTrack::Hive::PipeSeed::BlueprintDnasePeakCallSeed;

use strict;
use warnings;
use base ('ReseqTrack::Hive::PipeSeed::BasePipeSeed');
use ReseqTrack::Tools::Exception qw(throw warning);
use autodie;
use Data::Dump qw(dump);


sub create_seed_params {
  my ($self) = @_;
  my $options = $self->options;
  
  my $require_experiment_type = $options->{'require_experiment_type'} || undef;    
  
  throw('this module will only accept pipelines that work on the collection table')
      if $self->table_name ne 'collection';
      
  my $db = $self->db();
  my $sa = $db->get_SampleAdaptor;
  my $ea = $db->get_ExperimentAdaptor;
  my $ra = $db->get_RunAdaptor;
  my $ca = $db->get_CollectionAdaptor;
  my $fa = $db->get_FileAdaptor;


  $self->SUPER::create_seed_params();
  my @new_seed_params;
  
  SEED:
  foreach my $seed_params (@{$self->seed_params}) {
      my ($collection, $output_hash) = @$seed_params;

      my $collection_attributes = $collection->attributes;

      foreach my $collection_attribute (@$collection_attributes) {
        $output_hash->{ $collection_attribute->attribute_name } = $collection_attribute->attribute_value;    ## adding all collection attributes to pipeline
      }
     
      my $experiment_name = $collection->name;
      my $experiment = $ea->fetch_by_source_id( $experiment_name );
      

      $output_hash->{'experiment_source_id'} = $experiment_name;

   
        
      my $runs = $ra->fetch_by_experiment_id($experiment->dbID);  
      my $sample_id = $$runs[0]->sample_id;                       ## assuming all runs of an experiment are from same sample
      my $sample = $sa->fetch_by_dbID( $sample_id );
      $output_hash->{'sample_alias'} = $sample->sample_alias;     ## improve, make it generic
          
     push ( @new_seed_params, $seed_params );
  }  
dump( @new_seed_params );

  $self->seed_params(\@new_seed_params);  ## updating the seed param
  $db->dbc->disconnect_when_inactive(1);
}

1;


