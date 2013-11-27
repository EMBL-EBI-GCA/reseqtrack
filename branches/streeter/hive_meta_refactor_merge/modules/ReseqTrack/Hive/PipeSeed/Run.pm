
package ReseqTrack::Hive::PipeSeed::Run;

use strict;
use warnings;
require ReseqTrack::Hive::PipeSeed::Default;

sub output_params {
  return [
    @{ ReseqTrack::Hive::PipeSeed::Default::output_params()},
    qw(
    output_experiment_columns
    output_experiment_attributes
    output_sample_columns
    output_sample_attributes
    output_study_columns
    output_study_attributes
  )];
}


sub create_seed_params {
  my ($pipeline, $select_options, $output_params) = @_;

  my $output_experiment_columns = ref($output_params->{'output_experiment_columns'}) eq 'ARRAY' ? $output_params->{'output_experiment_columns'}
                      : defined $output_params->{'output_experiment_columns'} ? [$output_params->{'output_experiment_columns'}]
                      : [];
  my $output_experiment_attributes = ref($output_params->{'output_experiment_attributes'}) eq 'ARRAY' ? $output_params->{'output_experiment_attributes'}
                      : defined $output_params->{'output_experiment_attributes'} ? [$output_params->{'output_experiment_attributes'}]
                      : [];
  my $output_sample_columns = ref($output_params->{'output_sample_columns'}) eq 'ARRAY' ? $output_params->{'output_sample_columns'}
                      : defined $output_params->{'output_sample_columns'} ? [$output_params->{'output_sample_columns'}]
                      : [];
  my $output_sample_attributes = ref($output_params->{'output_sample_attributes'}) eq 'ARRAY' ? $output_params->{'output_sample_attributes'}
                      : defined $output_params->{'output_sample_attributes'} ? [$output_params->{'output_sample_attributes'}]
                      : [];
  my $output_study_columns = ref($output_params->{'output_study_columns'}) eq 'ARRAY' ? $output_params->{'output_study_columns'}
                      : defined $output_params->{'output_study_columns'} ? [$output_params->{'output_study_columns'}]
                      : [];
  my $output_study_attributes = ref($output_params->{'output_study_attributes'}) eq 'ARRAY' ? $output_params->{'output_study_attributes'}
                      : defined $output_params->{'output_study_attributes'} ? [$output_params->{'output_study_attributes'}]
                      : [];

  throw('this module will only accept pipelines that work on the run table')
      if $pipeline->table_name ne 'run';

  my $db = $pipeline->adaptor->db;

  my $sa = $db->get_SampleAdaptor;
  my $ea = $db->get_ExperimentAdaptor;
  my $sta = $db->get_StudyAdaptor;

  my $seed_params = ReseqTrack::Hive::PipeSeed::Default::create_seed_params($pipeline, $select_options, $output_params);
  foreach my $seed_params (@$seed_params) {
    my ($run, $output_hash) = @$seed_params;
    if (scalar @$output_sample_columns || scalar @$output_sample_attributes) {
      my $sample = $sa->fetch_by_dbID($run->sample_id);
      throw('did not get a sample with id '.$run->sample_id) if !$sample;
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
  return $seed_params;
}

1;
