package ReseqTrack::Tools::AttributeUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::Attribute;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use vars qw (@ISA  @EXPORT);

@ISA    = qw(Exporter);
@EXPORT = qw(create_attribute_for_object remove_outdated_attributes);

sub create_attribute_for_object {
	my ( $object, $name, $value, $units ) = @_;
	throw(
		    "Must pass create_attribute_for_object a ReseqTrack::HasHistory object "
			. "not "
			. $object )
		unless ( $object->isa("ReseqTrack::HasHistory") );
		
	throw("Can't create an attribute without a name "
			. $name
			. " or value "
			. $value . " for "
			. $object->name )
		unless ( $name && defined($value) && $value ne '' );
		
	my $attribute = ReseqTrack::Attribute->new(
		-table_name      => $object->object_table_name,
		-other_id        => $object->dbID,
		-attribute_name  => $name,
		-attribute_value => $value,
    -attribute_units => $units,
	);
	return $attribute;
}

sub remove_outdated_attributes {
	my ( $new_object, $old_object, $attribute_adaptor) =	@_;

	my %new_attr_names;
	map { $new_attr_names{ $_->attribute_name() } = 1 }
		@{ $new_object->attributes() };

	for my $attr ( @{ $old_object->attributes() } ) {
		my $key = $attr->attribute_name();
		if ( !exists $new_attr_names{$key} ) {
			$attribute_adaptor->delete($attr);
		}
	}
}

1;
