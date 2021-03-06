
package ReseqTrack::Hive::PipeSeed::BasePipeSeed;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);

sub sql_existing {
  my ($self)      = @_;
  my $table_name  = $self->table_name;
  my $dbID_name   = $table_name . '_id';
  my $pipeline_id = $self->pipeline_id;
  my $sql_existing =
      "SELECT pipeline_seed.seed_id FROM $table_name, pipeline_seed, hive_db"
    . " WHERE pipeline_seed.hive_db_id = hive_db.hive_db_id"
    . " AND hive_db.pipeline_id = $pipeline_id"
    . " AND $table_name.$dbID_name = pipeline_seed.seed_id"
    . " AND (pipeline_seed.is_running = 1"
    . " OR pipeline_seed.is_complete = 1"
    . " OR pipeline_seed.is_futile = 1)";
  return $sql_existing;
}

sub new {
  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;

  my ( $options, $pipeline ) = rearrange( [qw( OPTIONS PIPELINE )], @args );

  $self->options($options);
  $self->pipeline($pipeline);

  return $self;
}

sub pipeline {
  my ( $self, $pipeline ) = @_;
  if ($pipeline) {
    throw("pipeline is not a ReseqTrack::Pipeline object")
      if ref($pipeline) ne 'ReseqTrack::Pipeline';
    $self->{'pipeline'} = $pipeline;
  }
  return $self->{'pipeline'};
}

sub table_name {
  my ($self) = @_;
  return $self->pipeline->table_name;
}

sub pipeline_id {
  my ($self) = @_;
  return $self->pipeline->dbID;
}

sub db {
  my ($self) = @_;
  return $self->pipeline->adaptor->db;
}

sub dbID_name {
  my ($self) = @_;
  return $self->table_name . '_id';
}

sub options {
  my ( $self, $options ) = @_;
  if ($options) {
    throw("options is not a hashref") if ref($options) ne 'HASH';
    $self->{'options'} = $options;
  }
  return $self->{'options'} // {};
}

sub seed_params {
  my ( $self, $seed_params ) = @_;
  if ($seed_params) {
    $self->{'seed_params'} = $seed_params;
  }
  return $self->{'seed_params'} // [];
}

sub seed_sql_and_bind_variables {
  my ($self) = @_;

  my $table_name = $self->table_name();
  my $dbID_name  = $self->dbID_name();
  my $adaptor    = $self->db->get_adaptor_for_table( $self->table_name );

  my $options            = $self->options;
  my $require_columns    = $options->{'require_columns'} || {};
  my $require_attributes = $options->{'require_attributes'} || {};
  my $exclude_columns    = $options->{'exclude_columns'} || {};
  my $exclude_attributes = $options->{'exclude_attributes'} || {};

  my $sql_existing = $self->sql_existing();
  my $sql =
      "SELECT "
    . $adaptor->columns
    . " FROM $table_name "
    . " WHERE $dbID_name NOT IN ($sql_existing)";

  my @bind_values;
  while ( my ( $column, $values ) = each %$require_columns ) {
    if ( ref($values) eq 'ARRAY' ) {
      $sql .= " AND $table_name.$column IN ("
        . join( ',', map { '?' } @$values ) . ')';
      push( @bind_values, @$values );
    }
    else {
      $sql .= " AND $table_name.$column = ?";
      push( @bind_values, $values );
    }
  }
  while ( my ( $column, $values ) = each %$exclude_columns ) {
    if ( ref($values) eq 'ARRAY' ) {
      $sql .= " AND $table_name.$column NOT IN ("
        . join( ',', map { '?' } @$values ) . ')';
      push( @bind_values, @$values );
    }
    else {
      $sql .= " AND $table_name.$column != ?";
      push( @bind_values, $values );
    }
  }

  return ( $sql, @bind_values );
}

sub option_array {
  my ( $self, $option_name ) = @_;

  my $options = $self->options;

  return
      ref( $options->{$option_name} ) eq 'ARRAY' ? $options->{$option_name}
    : defined $options->{$option_name} ? [ $options->{$option_name} ]
    :                                    [];
}

sub create_seed_params {
  my ($self)            = @_;
  my $options           = $self->options;
  my $output_columns    = $self->option_array('output_columns');
  my $output_attributes = $self->option_array('output_attributes');

  my $table_name = $self->table_name();
  my $dbID_name  = $self->dbID_name();
  my $adaptor    = $self->db->get_adaptor_for_table( $self->table_name );

  my $require_columns    = $options->{'require_columns'}    || {};
  my $require_attributes = $options->{'require_attributes'} || {};
  my $exclude_columns    = $options->{'exclude_columns'}    || {};
  my $exclude_attributes = $options->{'exclude_attributes'} || {};

  my $db = $self->db();
  my ( $sql, @bind_values ) = $self->seed_sql_and_bind_variables();

  my $sth = $db->dbc->prepare($sql)
    or throw( "could not prepare $sql: " . $db->dbc->errstr );
  $sth->execute(@bind_values) or die "could not execute $sql: " . $sth->errstr;

  my @seed_params;

  my @remote_host_ids =
    map { $_->dbID } @{ $db->get_HostAdaptor->fetch_all_remote() };

SEED:
  while ( my $rowHashref = $sth->fetchrow_hashref ) {
    my $seed = $adaptor->object_from_hashref($rowHashref);

    # Do not run a pipeline on foreign files:
    next SEED
      if $table_name eq 'file' && grep { $seed->host_id == $_ }
      @remote_host_ids;
    if ( $table_name eq 'collection' && $seed->table_name eq 'file' ) {
      foreach my $file ( @{ $seed->others } ) {
        next SEED if grep { $file->host_id == $_ } @remote_host_ids;
      }
    }

    foreach my $attr_name ( keys %$require_attributes ) {
      my ($attribute) =
        grep { $_->attribute_name eq $attr_name } @{ $seed->attributes };
      next SEED if !$attribute;
      my $values = $require_attributes->{$attr_name};
      $values = ref($values) eq 'ARRAY' ? $values : [$values];
      next SEED if !grep { $attribute->attribute_value eq $_ } @$values;
    }
  ATTR:
    while ( my ( $attr_name, $values ) = each %$exclude_attributes ) {
      my ($attribute) =
        grep { $_->attribute_name eq $attr_name } @{ $seed->attributes };
      next ATTR if !$attribute;
      $values = ref($values) eq 'ARRAY' ? $values : [$values];
      next SEED if grep { $attribute->attribute_value eq $_ } @$values;
    }

    my %output_id;
    foreach my $column_name (@$output_columns) {
      throw("$column_name is not a valid column name for $table_name")
        if !exists( $rowHashref->{$column_name} );
      $output_id{$column_name} = $rowHashref->{$column_name};
    }
    if (@$output_attributes) {
      my $seed_attributes = $seed->attributes;
    ATTRIBUTE:
      foreach my $attribute_name (@$output_attributes) {
        my ($attribute) =
          grep { $_->attribute_name eq $attribute_name } @$seed_attributes;
        next ATTRIBUTE if !$attribute;
        $output_id{$attribute_name} = $attribute->attribute_value;
      }
    }
    push( @seed_params, [ $seed, \%output_id ] );
  }
  $self->seed_params( \@seed_params );
}

1;
