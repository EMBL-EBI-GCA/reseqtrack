package ReseqTrack::Hive::PipeSeed::CollectionFiles;

use strict;
use warnings;
use base ('ReseqTrack::Hive::PipeSeed::BasePipeSeed');

sub create_seed_params {
  my ($self) = @_;

  throw(
    'this module will only accept pipelines that work on the collection table')
    if $self->table_name ne 'collection';

  my $output_file_columns    = $self->option_array('output_file_columns');
  my $output_file_attributes = $self->option_array('output_file_attributes');

  my $db = $self->db();
  my $fa = $db->get_FileAdaptor;

  $self->SUPER::create_seed_params();

  for my $seed_params ( @{ $self->seed_params } ) {
    my ( $collection, $output_hash ) = @$seed_params;

    throw(
      'this module will only accept collections that work on file collections')
      if $collection->table_name ne 'file';

    my @files;
    $output_hash->{collection_files} = \@files;

    for my $file ( @{ $collection->others } ) {
      my $file_output_hash = {};
      push @files, $file_output_hash;

      foreach my $column_name (@$output_file_columns) {
        my $value;
        if ($column_name eq 'file_id') {
          $value = $file->dbID;
        }
        if ($column_name eq 'name') {
          $value = $file->name;
        }
        if ($column_name eq 'type') {
          $value = $file->type;
        }
        if ($column_name eq 'size') {
          $value = $file->size;
        }
        if ($column_name eq 'md5') {
          $value = $file->md5;
        }
        if ($column_name eq 'host_id') {
          $value = $file->host_id;
        }
        if ($column_name eq 'withdrawn') {
          $value = $file->withdrawn;
        }
        if ($column_name eq 'created') {
          $value = $file->created;
        }
        if ($column_name eq 'updated') {
          $value = $file->updated;
        }
        $file_output_hash->{$column_name} = $value;
      }

      if (@$output_file_attributes) {
        my $file_attributes = $file->attributes;
      ATTRIBUTE:
        foreach my $attribute_name (@$output_file_attributes) {
          my ($attribute) =
            grep { $_->attribute_name eq $attribute_name } @$file_attributes;
          next ATTRIBUTE if !$attribute;
          $output_hash->{$attribute_name} = $attribute->attribute_value;
        }
      }
    }

  }
}

1;
