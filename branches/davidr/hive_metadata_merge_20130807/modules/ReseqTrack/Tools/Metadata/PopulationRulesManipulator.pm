package ReseqTrack::Tools::Metadata::PopulationRulesManipulator;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw);
use base qw(ReseqTrack::Tools::Metadata::BaseMetadataManipulator);
use ReseqTrack::Tools::AttributeUtils
  qw(remove_outdated_attributes create_attribute_for_object);

use Data::Dumper;

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my $db               = $self->reseq_db();
  my $population_rules = $db->get_PopulationRuleAdaptor->fetch_all_in_order();

  $self->population_rules($population_rules);

  return $self;
}


sub convert_population {
  my ( $self, $sample ) = @_;
  my $attributes = $sample->attributes_hash;

  my $population = $attributes->{POPULATION}->attribute_value
    if ( $attributes->{POPULATION} );

  if ( !$population ) {
    my $adaptor = $self->era_db()->get_StudyAdaptor;
    my $studies = $adaptor->fetch_by_sample_id( $sample->source_id );
    my @names   = map { $_->title } @$studies;
    $population = join( ';', @names );
  }

  foreach my $pr ( @{ $self->population_rules } ) {
    if ( $pr->test_description($population) ) {
      $population = $pr->population;
      last;
    }
  }

  my $attr = create_attribute_for_object( $sample, 'POPULATION', $population );
  $sample->statistics( [$attr] );
}

sub manipulate_sample {
  my ( $self, $sample, $current_copy ) = @_;

  $self->convert_population($sample);
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
