package ReseqTrack::DBSQL::RunAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::LazyAdaptor);
use ReseqTrack::Run;

sub new {
	my ( $class, $db ) = @_;

	my $self = $class->SUPER::new($db);

	return $self;
}

sub column_mappings {
	my ( $self, $r ) = @_;

	throw("must be passed an object") unless ($r);

	return {
		run_id              => sub { $r->dbID(@_) },
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

sub fetch_by_ena_run_id {
	my ( $self, $ena_run_id ) = @_;
	my $runs = $self->fetch_by_column_name( 'ena_run_id', $ena_run_id );
	if (@$runs) {
		return @$runs;
	}
	return undef;
}

sub store {
	my ( $self, $run, $update ) = @_;
	my $existing_record = $self->fetch_by_ena_run_id( $run->ena_run_id );

	if ( $existing_record && $update ) {
		$run->dbID( $existing_record->dbID );
		return $self->update($run);
	}

	if ( !$existing_record ) {
		$self->SUPER::store( $run, $update );
	}
}

1;
