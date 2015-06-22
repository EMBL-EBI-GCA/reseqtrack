
package ReseqTrack::Hive::PipeSeed::BlueprintSeed;

use strict;
use warnings;
use base ('ReseqTrack::Hive::PipeSeed::BasePipeSeed');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);

sub sql_existing {
  my ($self) = @_;
  my $table_name = $self->table_name;
  my $dbID_name = $table_name . '_id';
  my $pipeline_id = $self->pipeline_id;
  my $sql_existing =
        "SELECT pipeline_seed.seed_id FROM $table_name, pipeline_seed, hive_db"
      . " WHERE pipeline_seed.hive_db_id = hive_db.hive_db_id"
      . " AND hive_db.pipeline_id = $pipeline_id"
      . " AND $table_name.$dbID_name = pipeline_seed.seed_id"
      . " AND (pipeline_seed.is_running = 1"
      .      " OR pipeline_seed.is_complete = 1"
      .      " OR pipeline_seed.is_futile = 1)";
  return $sql_existing;
}


sub create_seed_params {
  my ($self) = @_;
  my $options = $self->options;
  my $output_columns = ref($options->{'output_columns'}) eq 'ARRAY' ? $options->{'output_columns'}
                      : defined $options->{'output_columns'} ? [$options->{'output_columns'}]
                      : [];
  my $output_attributes = ref($options->{'output_attributes'}) eq 'ARRAY' ? $options->{'output_attributes'}
                      : defined $options->{'output_attributes'} ? [$options->{'output_attributes'}]
                      : [];

  my $require_study_columns = $options->{'require_study_columns'} || {};
  my $require_sample_columns = $options->{'require_sample_columns'} || {};
  my $require_experiment_columns = $options->{'require_experiment_columns'} || {};
  my $require_run_columns = $options->{'require_run_columns'} || {};
  my $require_study_attributes = $options->{'require_study_attributes'} || {};
  my $require_sample_attributes = $options->{'require_sample_attributes'} || {};
  my $require_experiment_attributes = $options->{'require_experiment_attributes'} || {};
  my $require_run_attributes = $options->{'require_run_attributes'} || {};
  my $exclude_study_columns = $options->{'exclude_study_columns'} || {};
  my $exclude_sample_columns = $options->{'exclude_sample_columns'} || {};
  my $exclude_experiment_columns = $options->{'exclude_experiment_columns'} || {};
  my $exclude_run_columns = $options->{'exclude_run_columns'} || {};
  my $exclude_study_attributes = $options->{'exclude_study_attributes'} || {};
  my $exclude_sample_attributes = $options->{'exclude_sample_attributes'} || {};
  my $exclude_experiment_attributes = $options->{'exclude_experiment_attributes'} || {};
  my $exclude_run_attributes = $options->{'exclude_run_attributes'} || {};

  my $db = $self->db();
  my $table_name = $self->table_name();
  my $dbID_name = $self->dbID_name();
  my $adaptor = $self->db->get_adaptor_for_table($self->table_name);

  my $sql_existing = $self->sql_existing();
  
  my @all_columns;

  my $study_adaptor = $db->get_StudyIDAdaptor;
  my $study_columns = $study_adaptor->columns;
  $study_columns = map{ "study.".$_ ) @{$study_columns};
  push @all_columns, $study_columns;

  my $sample_adaptor = $db->get_SampleAdaptor;
  my $sample_columns = $sample_adaptor->columns;
  $sample_columns = map{ "sample.".$_ ) @{$sample_columns};
  push @all_columns, $sample_columns;

  my $experiment_adaptor = $db->get_ExperimentAdaptor;
  my $experiment_columns = $experiment_adaptor->columns;
  $experiment_columns = map{ "experiment.".$_ ) @{$experiment_columns};
  push @all_columns, $experiment_columns;

  my $run_adaptor = $db->get_RunAdaptor;
  my $run_columns = $run_adaptor->columns;
  $run_columns = map{ "run.".$_ ) @{$run_columns}; 
  push @all_columns, $run_columns;
  

 
  my $sql = "SELECT $all_columns FROM study, sample, experiment, run ".
          . "WHERE run."   
     ;
  
  my $sql = "SELECT ".$adaptor->columns." FROM $table_name "
        . " WHERE $dbID_name NOT IN ($sql_existing)";
  my @bind_values;
  while (my ($column, $values) = each %$require_columns) {
    if (ref($values) eq 'ARRAY') {
      $sql .= " AND $table_name.$column IN (" . join(',', map {'?'} @$values) . ')';
      push(@bind_values, @$values);
    }
    else {
      $sql .= " AND $table_name.$column = ?";
      push(@bind_values, $values);
    }
  }
  while (my ($column, $values) = each %$exclude_columns) {
    if (ref($values) eq 'ARRAY') {
      $sql .= " AND $table_name.$column NOT IN (" . join(',', map {'?'} @$values) . ')';
      push(@bind_values, @$values);
    }
    else {
      $sql .= " AND $table_name.$column != ?";
      push(@bind_values, $values);
    }
  }
  my $sth = $db->dbc->prepare($sql) or throw("could not prepare $sql: ".$db->dbc->errstr);
  $sth->execute(@bind_values) or die "could not execute $sql: ".$sth->errstr;

  my @seed_params;

  my @remote_host_ids = map {$_->dbID} @{$db->get_HostAdaptor->fetch_all_remote()};

  SEED:
  while (my $rowHashref = $sth->fetchrow_hashref) {
    my $seed = $adaptor->object_from_hashref($rowHashref);

    # Do not run a pipeline on foreign files:
    next SEED if $table_name eq 'file' && grep {$seed->host_id == $_} @remote_host_ids;
    if ($table_name eq 'collection' && $seed->table_name eq 'file') {
      foreach my $file (@{$seed->others}) {
        next SEED if grep {$file->host_id == $_} @remote_host_ids;
      }
    }

    foreach my $attr_name (keys %$require_attributes) {
      my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$seed->attributes};
      next SEED if !$attribute;
      my $values = $require_attributes->{$attr_name};
      $values = ref($values) eq 'ARRAY' ? $values : [$values];
      next SEED if !grep {$attribute->attribute_value eq $_} @$values;
    }
    ATTR:
    while (my ($attr_name, $values) = each %$exclude_attributes) {
      my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$seed->attributes};
      next ATTR if !$attribute;
      $values = ref($values) eq 'ARRAY' ? $values : [$values];
      next SEED if grep {$attribute->attribute_value eq $_} @$values;
    }

    my %output_id;
    foreach my $column_name (@$output_columns) {
      throw("$column_name is not a valid column name for $table_name") if !exists($rowHashref->{$column_name});
      $output_id{$column_name} = $rowHashref->{$column_name};
    }
    if (@$output_attributes) {
      my $seed_attributes = $seed->attributes;
      ATTRIBUTE:
      foreach my $attribute_name (@$output_attributes) {
        my ($attribute) = grep {$_->attribute_name eq $attribute_name} @$seed_attributes;
        next ATTRIBUTE if !$attribute;
        $output_id{$attribute_name} = $attribute->attribute_value;
      }
    }
    push(@seed_params, [$seed, \%output_id]);
  }
  $self->seed_params(\@seed_params);
}

1;
