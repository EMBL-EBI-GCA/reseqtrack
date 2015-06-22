package ReseqTrack::Tools::Metadata::G1KManipulator;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw);
use base qw(ReseqTrack::Tools::Metadata::BaseMetadataManipulator);
use ReseqTrack::Tools::AttributeUtils
  qw(remove_outdated_attributes create_attribute_for_object);
  
sub validate_sample_name {
  my ($self,$sample) = @_;
  if (! $sample->sample_alias =~ m/^[A-Z]{2}\d{5}$/){
    throw("Sample alias does not match expected format ".$sample->source_id." ".$sample->sample_alias);
  }
}

sub update_sample_name {
  my ($self,$sample) = @_;
  
  my $sample_name = $sample->sample_alias;
  $sample_name =~ s/^GM/NA/;
  $sample_name = uc($sample_name);
  
  $sample->sample_alias($sample_name);
}

sub manipulate_sample {
  my ( $self, $sample, $current_copy ) = @_;
  $self->update_sample_name($sample);
  $self->validate_sample_name($sample);
}

sub population_rules {
  my ( $self, $rules ) = @_;
  if ($rules) {
    throw("rules must be an arrayref") if ( ref $rules ne 'ARRAY' );
    $self->{population_rules} = $rules;
  }
  return $self->{population_rules};
}
1;
