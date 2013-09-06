
package ReseqTrack::Hive::PipeSeed::Default;

use strict;
use warnings;

sub output_params {
  return [qw(
    output_columns
    output_attributes
  )];
}

sub sql_existing {
  my ($pipeline) = @_;
  my $table_name = $pipeline->table_name;
  my $dbID_name = $table_name . '_id';
  my $pipeline_id = $pipeline->dbID;
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


#select_options affect what gets selected as a pipeline seed.  These are defined in the ReseqTrack pipeline table
#output_params affect what parameters get written to the output id. Defined by the pipeline configuration.

sub create_seed_params {
  my ($pipeline, $select_options, $output_params) = @_;
  my $output_columns = ref($output_params->{'output_columns'}) eq 'ARRAY' ? $output_params->{'output_columns'}
                      : defined $output_params->{'output_columns'} ? [$output_params->{'output_columns'}]
                      : [];
  my $output_attributes = ref($output_params->{'output_attributes'}) eq 'ARRAY' ? $output_params->{'output_attributes'}
                      : defined $output_params->{'output_attributes'} ? [$output_params->{'output_attributes'}]
                      : [];

  my $require_columns = $select_options->{'require_columns'} || {};
  my $require_attributes = $select_options->{'require_attributes'} || {};
  my $exclude_columns = $select_options->{'exclude_columns'} || {};
  my $exclude_attributes = $select_options->{'exclude_attributes'} || {};

  my $db = $pipeline->adaptor->db;

  my $table_name = $pipeline->table_name;
  my $dbID_name = $table_name . '_id';

  my $adaptor = $db->get_adaptor_for_table($table_name);
  my $sql_existing = sql_existing($pipeline);
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

  SEED:
  while (my $rowHashref = $sth->fetchrow_hashref) {
    my $seed = $adaptor->object_from_hashref($rowHashref);
    while (my ($attr_name, $values) = each %$require_attributes) {
      my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$seed->attributes};
      next SEED if !$attribute;
      $values = ref($values) eq 'ARRAY' ? $values : [$values];
      next SEED if !grep {$attribute->value eq $_} @$values;
    }
    ATTR:
    while (my ($attr_name, $values) = each %$exclude_attributes) {
      my ($attribute) = grep {$_->attribute_name eq $attr_name} @{$seed->attributes};
      next ATTR if !$attribute;
      $values = ref($values) eq 'ARRAY' ? $values : [$values];
      next SEED if grep {$attribute->value eq $_} @$values;
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
  return \@seed_params;
}

1;
