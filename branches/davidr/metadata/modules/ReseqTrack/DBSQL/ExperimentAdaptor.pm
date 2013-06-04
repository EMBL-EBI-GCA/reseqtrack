package ReseqTrack::DBSQL::ExperimentAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::LazyAdaptor);

use ReseqTrack::Experiment;


sub new {
	my ( $class, $db ) = @_;

	my $self = $class->SUPER::new($db);

	return $self;
}

sub column_mappings {
	my ( $self, $e ) = @_;

	throw("must be passed an object") unless ($e);

	return {
		experiment_id         => sub { $e->dbID(@_) },
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
		source_id             => sub { $e->source_id(@_) },
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
	return pop @{ $self->fetch_by_column_name( "source_id", $source_id ) };
}

sub store {
	my ( $self, $experiment, $update ) = @_;
	my $existing_record = $self->fetch_by_dbID( $experiment->dbID ) if ($experiment->dbID);

	if ( $existing_record && $update ) {
		$experiment->dbID( $existing_record->dbID );
		return $self->update($experiment);
	}

	if ( !$existing_record ) {
		$self->SUPER::store( $experiment, $update );
	}
}

1;
