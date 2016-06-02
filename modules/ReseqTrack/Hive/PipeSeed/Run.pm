
package ReseqTrack::Hive::PipeSeed::Run;

use strict;
use warnings;
use base ('ReseqTrack::Hive::PipeSeed::BasePipeSeed');


sub create_seed_params {
  my ($self) = @_;
  my $options = $self->options;

  my $output_experiment_columns = $self->option_array('output_experiment_columns');
  my $output_experiment_attributes = $self->option_array('output_experiment_attributes');
  my $output_sample_columns = $self->option_array('output_sample_columns');
  my $output_sample_attributes = $self->option_array('output_sample_attributes');
  my $output_study_columns = $self->option_array('$output_study_columns');
  my $output_study_attributes = $self->option_array('$$output_study_attributes');

  throw('this module will only accept pipelines that work on the run table')
      if $self->table_name ne 'run';

  my $db = $self->db();
  my $sa = $db->get_SampleAdaptor;
  my $ea = $db->get_ExperimentAdaptor;
  my $sta = $db->get_StudyAdaptor;

  $self->SUPER::create_seed_params();

  foreach my $seed_params (@{$self->seed_params}) {
    my ($run, $output_hash) = @$seed_params;
    if (scalar @$output_sample_columns || scalar @$output_sample_attributes) {
      my $sample = $run->experiment->sample;
      throw('did not get a sample for run '.$run->run_source_id) if !$sample;
      foreach my $column_name (@$output_sample_columns) {
        $output_hash->{$column_name} = &{$sa->column_mappings($sample)->{$column_name}}();
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
        $output_hash->{$column_name} = &{$ea->column_mappings($experiment)->{$column_name}}();
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
          $output_hash->{$column_name} = &{$sta->column_mappings($study)->{$column_name}}();
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
}

1;
