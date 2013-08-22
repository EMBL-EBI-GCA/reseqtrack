
package ReseqTrack::Hive::PipeSeed::Default;

use strict;
use warnings;

sub create_seed_params {
  my ($seed_factory, $pipeline, $select_options) = @_;
  my $output_columns = $seed_factory->param_to_flat_array('output_columns');
  my $output_attributes = $seed_factory->param_to_flat_array('output_attributes');

  my $require_columns = $select_options->{'require_columns'} || {};
  my $require_attributes = $select_options->{'require_attributes'} || {};
  my $exclude_columns = $select_options->{'exclude_columns'} || {};
  my $exclude_attributes = $select_options->{'exclude_attributes'} || {};

  my $db = $pipeline->adaptor->db;

  my $table_name = $pipeline->table_name;
  my $dbID_name = $table_name . '_id';

  my $adaptor = $db->get_adaptor_for_table($table_name);
  my $sql_existing = $seed_factory->sql_existing($pipeline);
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
    foreach my $output_column (@$output_columns) {
      my ($output_name, $column_name) = ref($output_column) eq 'HASH' ? %$output_column : ($output_column, $output_column);
      throw("$column_name is not a valid column name for $table_name") if !exists($rowHashref->{$column_name});
      $output_id{$output_name} = $rowHashref->{$column_name};
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
