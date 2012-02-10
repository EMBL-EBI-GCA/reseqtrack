package ReseqTrack::Tools::StatisticsUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::Statistic;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(create_statistic_for_object);


sub create_statistic_for_object{
  my ($object, $name, $value) = @_;
  throw("Must pass create_statistic_for_object a ReseqTrack::HasHistory object ".
        "not ".$object) unless($object->isa("ReseqTrack::HasHistory"));
  throw("Can't create a statistic without a name ".$name." or value ".$value.
	" for ".$object->name) unless($name && defined($value) && $value ne '');
  my $statistic = ReseqTrack::Statistic->new(
    -table_name => $object->object_table_name,
    -other_id => $object->dbID,
    -attribute_name => $name,
    -attribute_value => $value
      );
  return $statistic;
}

1;
