package ReseqTrack::DBSQL::ERA::BaseEraAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::LazyAdaptor);
use ReseqTrack::Tools::Exception qw(throw);

sub store {
	throw("We do not update ERAPRO");
}

sub update {
	throw("We do not update ERAPRO");
}

sub column_mappings {
	 throw("Need to implement the column_mappings method");
}

sub object_class {
	throw("Need to specify the class this adaptor returns ");
}

sub xml_column {
	throw("Need to implement xml_column");
}

sub handle_xml {
	throw("Need to implement handle_xml");
}

sub columns {
	my ($self) = @_;
	my @columns = $self->SUPER::columns();
	my $xml_column = $self->xml_column();
	my $table_name = $self->table_name();
	push @columns, "$table_name.$xml_column" if ($xml_column);
	return @columns;
}

sub object_from_hash_ref {
	my ( $self, $hashref ) = @_;
	
	my $obj = $self->SUPER::object_from_hash_ref($hashref);	
	
	my $xml_col = $self->xml_column();
	if ($xml_col){
		my $xml = $hashref->{$self->xml_col()};
		handle_xml($xml,$obj);
	}
	
	return $obj;
}

1;