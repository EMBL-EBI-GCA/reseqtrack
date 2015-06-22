package ReseqTrack::Tools::Metadata::EnaReadInfo;

use base qw(ReseqTrack::Tools::Metadata::BaseMetadataManipulator);
use ReseqTrack::Tools::AttributeUtils qw(create_attribute_for_object);
use Data::Dumper;

sub manipulate_run {
  my ( $self, $run, $current_copy ) = @_;
  my $run_adaptor = $self->era_db->get_RunAdaptor();

  my $changes_made = 0;
  my $attribute_hash= {};
  if ($current_copy)  {
    $attributes_hash = $current_copy->attributes_hash();
  }
  
  my $stats = $run_adaptor->get_run_stats( $run->source_id );

  $changes_made +=
    $self->add_attribute_if_changed( $run, $attributes_hash, 'READ_COUNT',
    $stats->{READ_COUNT} );
  $changes_made +=
    $self->add_attribute_if_changed( $run, $attributes_hash, 'BASE_COUNT',
    $stats->{BASE_COUNT} );
  $changes_made +=
    $self->add_attribute_if_changed( $run, $attributes_hash, 'FASTQ_DATE',
    $stats->{FASTQ_DATE} );
  $changes_made +=
    $self->add_attribute_if_changed( $run, $attributes_hash, 'FASTQ_ERROR',
    $stats->{FASTQ_ERROR} );

  return $changes_made;
}

sub add_attribute_if_changed {
  my ( $self, $run, $current_attributes, $name, $value ) = @_;
  
  return 0 if (! defined $value);

  my $update_required = 0;
  
  if ( $current_attributes->{$name} ) {
    $update_required = 1 if ( $current_attributes->{$name}->attribute_value() ne $value );
  }
  else {
    $update_required = 1;
  }

  if ($update_required) {
    my $attr = create_attribute_for_object( $run, $name, $value );
    $run->statistics( [$attr] );
  }
  
  return $update_required;
}

1;
