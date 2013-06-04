package ReseqTrack::DBSQL::ERA::RunAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::ERA::BaseEraAdaptor);
use ReseqTrack::Run;

sub column_mappings {
	my ( $self, $r ) = @_;

	throw("must be passed an object") unless ($r);

	return {
		run_id              => sub { $r->source_id(@_) },
		experiment_id       => sub { $r->experiment_id(@_) },
		sample_id           => sub { $r->sample_id(@_) },
		run_alias           => sub { $r->run_alias(@_) },
		status              => sub { $r->status(@_) },
		md5                 => sub { $r->md5(@_) },
		center_name         => sub { $r->center_name(@_) },
		run_center_name     => sub { $r->run_center_name(@_) },
		instrument_platform => sub { $r->instrument_platform(@_) },
		instrument_model    => sub { $r->instrument_model(@_) },
	};
}

sub object_class {
	return 'ReseqTrack::Run';
}

sub table_name {
	return "run";
}

sub fetch_by_source_id {
	my ( $self, $source_id ) = @_;
	return pop @{ $self->fetch_by_column_name( "run_id", $source_id ) };
}

sub xml_column {
	return "";
}

1;
