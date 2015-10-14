package ReseqTrack::Hive::PipeSeed::ChipPeakCallSeed;

use strict;
use warnings;
use base ('ReseqTrack::Hive::PipeSeed::BasePipeSeed');
use ReseqTrack::Tools::Exception qw(throw warning);
use autodie;
use Data::Dump qw(dump);

sub default_options {
  return {
    input_prefix => 'Input_file',
    exp_type_attribute_name => 'EXPERIMENT_TYPE',
    require_experiment_type =>  'Input',
  };
}

sub create_seed_params {
  my ($self) = @_;
  my $options = $self->options;

  my $metadata_file = $options->{'metadata_file'};
  throw('require metadata file') unless $metadata_file;
  
  my $path_names_array = $options->{'path_names_array'} ? $options->{'path_names_array'} : undef;
  my $metadata_hash    = _get_metadata_hash( $metadata_file, 'EXPERIMENT_ID' );

  my $require_experiment_type = $options->{'require_experiment_type'} || undef;    
  my $non_match_input = $options->{'non_match_input'} || undef;
  my $bam_collection_type = $options->{'bam_collection_type'} || undef;
  my $input_prefix = $options->{'input_prefix'} || undef;
  my $exp_type_attribute_name = $options->{'exp_type_attribute_name'} || undef;
  
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
  my $input_hash = _read_non_match_input( $non_match_input ) if $non_match_input;
  
  SEED:
  foreach my $seed_params (@{$self->seed_params}) {
      my ($collection, $output_hash) = @$seed_params;

      my $collection_attributes = $collection->attributes;

      my $experiment_name = $collection->name;
      throw("$experiment_name not present in $metadata_file") unless exists ( $$metadata_hash{$experiment_name} );
      my $metadata_path_hash   = _get_path_hash( $experiment_name, $metadata_hash, $path_names_array );

      foreach my $path_name ( keys %{$metadata_path_hash} ){
        my $path_value = $$metadata_path_hash{$path_name};
        $output_hash->{$path_name} = $path_value;
      }

      foreach my $collection_attribute (@$collection_attributes) {
        $output_hash->{ $collection_attribute->attribute_name } = $collection_attribute->attribute_value;    ## adding all collection attributes to pipeline
      }
      
      my $experiment = $ea->fetch_by_source_id( $experiment_name );
      my $experiment_attributes = $experiment->attributes;
      my ($attribute) = grep {$_->attribute_name eq $exp_type_attribute_name} @$experiment_attributes;
      next SEED if $attribute->attribute_value eq $require_experiment_type;          ## not seeding pipeline for input experiments
      
      my $attribute_name = $attribute->attribute_name;
      my $attribute_value = $attribute->attribute_value;

      $attribute_value=~ s/Histone\s+//g
         if( $attribute_name eq 'EXPERIMENT_TYPE' );    ## fix for blueprint ChIP file name

      $attribute_value=~ s/ChIP-Seq\s+//g
         if( $attribute_name eq 'EXPERIMENT_TYPE' );    ## fix for blueprint ChIP file name

    #  $attribute_value=~ s/Chromatin\sAccessibility/Dnase/
    #     if( $attribute_name eq 'EXPERIMENT_TYPE' );    ## fix for blueprint Dnase file name 

      $output_hash->{$attribute_name} = $attribute_value;
      $output_hash->{'experiment_source_id'} = $experiment_name;


      my $broad = assign_peak_call_type( $attribute_value )
                      if( $attribute_name eq 'EXPERIMENT_TYPE' );
   
      $output_hash->{'broad'} = $broad;    ## setting parameter for peak calling 
 
      throw("did not get an experiment wiith $experiment_name") if !$experiment;
        
      my $input_file;
      my $runs = $ra->fetch_by_experiment_id($experiment->dbID);  
      my $sample_id = $$runs[0]->sample_id;                    ## assuming all runs of an experiment are from same sample
      my $sample = $sa->fetch_by_dbID( $sample_id );
      $output_hash->{'sample_alias'} = $sample->sample_alias;
        
      if ( exists( $$input_hash{ $experiment_name } )) {
          $input_file = $$input_hash{ $experiment_name };
      }
      else {   
    #    my $runs = $ra->fetch_by_experiment_id($experiment->dbID);  
    #    my $sample_id = $$runs[0]->sample_id;                    ## assuming all runs of an experiment are from same sample
    #    my $sample = $sa->fetch_by_dbID( $sample_id );
    #    $output_hash->{'sample_alias'} = $sample->sample_alias;

        my $all_experiments = $ea->fetch_by_sample_id( $sample_id );

        EXP:
        foreach my $exp( @$all_experiments ){
          my $experiment_attributes = $exp->attributes;
          my ($attribute) = grep {$_->attribute_name eq $exp_type_attribute_name} @$experiment_attributes;
          next EXP unless $attribute->attribute_value eq $require_experiment_type;                             ## look into next experiment
            
          my $input_experiment_name =  $exp->source_id;

          throw('require bam_collection_type') if !$bam_collection_type;
          my $collection = $ca->fetch_by_name_and_type( $input_experiment_name, $bam_collection_type );
          next EXP if !$collection;

          throw("got multiple file for collection $input_experiment_name and type $bam_collection_type") if scalar @{$collection->other_ids} >1;

          my $file = $fa->fetch_by_dbID( $collection->other_ids->[0]);
          next EXP if !$file;

          throw("got multiple input file for experiment $experiment_name") if $input_file;  ## expecting single input file, unless mentioned in the non_match_input file
          $input_file = $file->name; 
          
        }
     }
     next SEED  if !$input_file;  ## input may not be available yet
     
     throw('no input file prefix') if !$input_prefix;
     $output_hash->{$input_prefix} =    $input_file;
     push ( @new_seed_params, $seed_params );
  }  
dump( @new_seed_params );

  $self->seed_params(\@new_seed_params);  ## updating the seed param
  $db->dbc->disconnect_when_inactive(1);
}
  
sub _read_non_match_input {
  my ($file) =@_;
  my %input_hash;
    
  open my $fh, '<', $file;
  while( <$fh> ){
    chomp;
    next if (/^#/);
    my ($experiment, $input_file_path) = split '\t';
    $input_hash{ $experiment } = $input_file_path;
  }
  return \%input_hash;
}

sub assign_peak_call_type {
  my ( $exp_type ) = @_;
  my $broad = 0;

  my %chip_marks_hash = ( 'H3K27me3' => 1,
                          'H3K36me3' => 1,
                          'H3K4me1'  => 1,
                          'H3K9me3'  => 1,
                        ); 
  
 $broad = 1 if ( exists ( $chip_marks_hash{ $exp_type }));
 return $broad;
}

=head1 _get_path_hash

Returns a metadata hash from metadat hash. Inputs are key metadata id, metadata hash and an array of parameters (optional).

=cut

sub _get_path_hash {
  my ( $key_id, $metadata_hash, $path_names_array ) = @_;
  my $path_hash;

  throw("$key_id not found in metadata file") unless exists $$metadata_hash{ $key_id };
 
  $$metadata_hash{ $key_id }{SAMPLE_DESC_1} = "NO_TISSUE" 
    if ( $$metadata_hash{ $key_id }{SAMPLE_DESC_1} eq "-" );
 
  $$metadata_hash{ $key_id }{SAMPLE_DESC_2} = "NO_SOURCE" 
    if ( $$metadata_hash{ $key_id }{SAMPLE_DESC_2} eq "-" );
  
  $$metadata_hash{ $key_id }{SAMPLE_DESC_3} = "NO_CELL_TYPE" 
    if ( $$metadata_hash{ $key_id }{SAMPLE_DESC_3} eq "-" );

  if ( scalar @$path_names_array >= 1 ){
    my @uc_path_names_array = map{ uc($_) } @$path_names_array;
    my $key_metadata_hash   = $$metadata_hash{ $key_id };

    foreach my $key( @uc_path_names_array ){
      throw("$key in not present in metadata") unless exists $$key_metadata_hash{ $key };
    }

    my @path_name_values    = @$key_metadata_hash{ @uc_path_names_array };
    @path_name_values       = map{ s/[\s=\/\\;()]/_/g; $_; }@path_name_values;
    @path_name_values       = map{ s/_+/_/g; $_; }@path_name_values;
    @path_name_values       = map{ s/_$//g; $_; }@path_name_values;

    @$path_hash{ @$path_names_array } = @path_name_values;
  }
  else {
    $path_hash = $$metadata_hash{ $key_id };
  }
  return $path_hash;
}


=head1 _get_metadata_hash

Returns metadata hash from an index file keyed by any selected field.

=cut

sub _get_metadata_hash {
my ( $file, $key_string ) = @_;
  open my $fh, '<', $file;
  my @header;
  my %data;
  my $key_index = undef;

  while ( <$fh> ) {
    chomp;
    next if m/^#/;
    my @vals = split "\t", $_;

    if ( @header ) {
      throw("$key_string not found in $file") unless $key_index >= 0;
      $data { $vals[$key_index] }{ $header[$_] } = $vals[$_] for 0..$#header;
    }
    else {
      @header = map { uc($_) } @vals;
      my @key_index_array = grep{ $header[$_] eq $key_string } 0..$#header;
      $key_index = $key_index_array[0];
    }
  }
  return \%data;
  close( $fh );
}


1;


