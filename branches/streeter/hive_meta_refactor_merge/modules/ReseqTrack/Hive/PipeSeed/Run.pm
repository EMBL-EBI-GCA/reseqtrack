
package ReseqTrack::Hive::PipeSeed::Run;

use strict;
use warnings;
require ReseqTrack::Hive::PipeSeed::Default;

sub create_seed_params {
  my ($seed_factory, $pipeline, $select_options) = @_;
  my $output_run_columns = $seed_factory->param_to_flat_array('output_run_columns');
  my $output_run_attributes = $seed_factory->param_to_flat_array('output_run_attributes');
  my $output_sample_columns = $seed_factory->param_to_flat_array('output_sample_columns');
  my $output_sample_attributes = $seed_factory->param_to_flat_array('output_sample_attributes');
  my $output_experiment_columns = $seed_factory->param_to_flat_array('output_experiment_columns');
  my $output_experiment_attributes = $seed_factory->param_to_flat_array('output_experiment_attributes');
  my $output_study_columns = $seed_factory->param_to_flat_array('output_study_columns');
  my $output_study_attributes = $seed_factory->param_to_flat_array('output_study_attributes');

  throw('this module will only accept pipelines that work on the run table')
      if $pipeline->table_name ne 'run';

  my $db = $pipeline->adaptor->db;

  my $sa = $db->get_SampleAdaptor;
  my $ea = $db->get_ExperimentAdaptor;
  my $sta = $db->get_StudyAdaptor;


  my $seed_params = ReseqTrack::Hive::PipeSeed::Default::create_seed_params($seed_factory, $pipeline, $select_options);
  foreach my $seed_params (@$seed_params) {
    my ($run, $output_hash) = @$seed_params;
    if (scalar @$output_sample_columns || scalar @$output_sample_attributes) {
      my $sample = $sa->fetch_by_dbID($run->sample_id);
      throw('did not get a sample with id '.$run->sample_id) if !$sample;
      foreach my $column_name (@$output_sample_columns) {
        $output_hash->{$column_name} = &{$sa->column_mappings($sample)->{$column_name}};
      }
      if (@$output_sample_attributes) {
        my $sample_attributes = $sample->attributes;
        ATTRIBUTE:
        foreach my $attribute_name (@$output_sample_attributes) {
          my ($attribute) = grep {$_->attribute_name eq $attribute_name} @$sample_attributes;
          next ATTRIBUTE if !$attribute;
          $output_hash->{$attribute_name} = $attribute->attribute_value;
        }
      }
    }

    if (scalar @$output_experiment_columns || scalar @$output_experiment_attributes
          || scalar @$output_study_columns || scalar @$output_study_attributes) {
      my $experiment = $ea->fetch_by_dbID($run->experiment_id);
      throw('did not get an experiment with id '.$run->experiment_id) if !$experiment;
      foreach my $column_name (@$output_experiment_columns) {
        $output_hash->{$column_name} = &{$sa->column_mappings($experiment)->{$column_name}};
      }
      if (@$output_experiment_attributes) {
        my $experiment_attributes = $experiment->attributes;
        ATTRIBUTE:
        foreach my $attribute_name (@$output_experiment_attributes) {
          my ($attribute) = grep {$_->attribute_name eq $attribute_name} @$experiment_attributes;
          next ATTRIBUTE if !$attribute;
          $output_hash->{$attribute_name} = $attribute->attribute_value;
        }
      }
      if (scalar @$output_study_columns || scalar @$output_study_attributes) {
        my $study = $sta->fetch_by_dbID($experiment->study_id);
        throw('did not get a study with id '.$experiment->study_id) if !$study;
        foreach my $column_name (@$output_study_columns) {
          $output_hash->{$column_name} = &{$sta->column_mappings($study)->{$column_name}};
        }
        if (@$output_study_attributes) {
          my $study_attributes = $study->attributes;
          ATTRIBUTE:
          foreach my $attribute_name (@$output_study_attributes) {
            my ($attribute) = grep {$_->attribute_name eq $attribute_name} @$study_attributes;
            next ATTRIBUTE if !$attribute;
            $output_hash->{$attribute_name} = $attribute->attribute_value;
          }
        }
      }
    }
  }
  return $seed_params;
}

1;
