package ReseqTrack::DBSQL::ERA::SampleAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::ERA::BaseEraAdaptor);
use ReseqTrack::Sample;

sub column_mappings {
	my ( $self, $s ) = @_;

	throw("must be passed an object") unless ($s);

	return {
		sample_id       => sub { $s->source_id(@_) },
		status          => sub { $s->status(@_) },
		md5             => sub { $s->md5(@_) },
		center_name     => sub { $s->center_name(@_) },
		sample_alias    => sub { $s->sample_alias(@_) },
		tax_id          => sub { $s->tax_id(@_) },
		scientific_name => sub { $s->scientific_name(@_) },
		common_name     => sub { $s->common_name(@_) },
		anonymized_name => sub { $s->anonymized_name(@_) },
		individual_name => sub { $s->individual_name(@_) },
		sample_title    => sub { $s->sample_title(@_) },
	};
}

sub object_class {
	return 'ReseqTrack::Sample';
}

sub table_name {
	return "sample";
}

sub fetch_by_source_id {
	my ( $self, $source_id ) = @_;
	return pop @{ $self->fetch_by_column_name( "sample_id", $source_id ) };
}

sub xml_column {
	return "sample_xml"
}

sub handle_xml {
	#TODO
}


1;
