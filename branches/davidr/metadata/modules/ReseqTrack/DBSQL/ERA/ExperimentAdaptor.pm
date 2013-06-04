package ReseqTrack::DBSQL::ERA::ExperimentAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::ERA::BaseEraAdaptor);

use ReseqTrack::Experiment;

sub column_mappings {
	my ( $self, $e ) = @_;

	throw("must be passed an object") unless ($e);

	return {
		experiment_id         => sub { $e->source_id(@_) },
		study_id              => sub { $e->study_id(@_) },
		status                => sub { $e->status(@_) },
		md5                   => sub { $e->md5(@_) },
		center_name           => sub { $e->center_name(@_) },
		experiment_alias      => sub { $e->experiment_alias(@_) },
		instrument_platform   => sub { $e->instrument_platform(@_) },
		instrument_model      => sub { $e->instrument_model(@_) },
		library_layout        => sub { $e->library_layout(@_) },
		library_name          => sub { $e->library_name(@_) },
		library_strategy      => sub { $e->library_strategy(@_) },
		library_source        => sub { $e->library_source(@_) },
		library_selection     => sub { $e->library_selection(@_) },
		paired_nominal_length => sub { $e->paired_nominal_length(@_) },
		paired_nominal_sdev   => sub { $e->paired_nominal_sdev(@_) },
	};
}

sub object_class {
	return 'ReseqTrack::Experiment';
}

sub table_name {
	return "experiment";
}

sub fetch_by_source_id {
	my ( $self, $source_id ) = @_;
	return pop @{ $self->fetch_by_column_name( "experiment_id", $source_id ) };
}

sub xml_column {
	return "experiment_xml";
}

sub handle_xml {
	#TODO
}

1;
