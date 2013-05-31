package ReseqTrack::DBSQL::SampleAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::LazyAdaptor);
use ReseqTrack::Sample;

sub new {
	my ( $class, $db ) = @_;

	my $self = $class->SUPER::new($db);

	return $self;
}

sub column_mappings {
	my ( $self, $s ) = @_;

	throw("must be passed an object") unless ($s);

	return {
		sample_id       => sub { $s->dbID(@_) },
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
	return 'ReseqTrack::Run';
}

sub table_name {
	return "run";
}

sub fetch_by_ena_sample_id {
	my ( $self, $ena_sample_id ) = @_;
	my $runs = $self->fetch_by_column_name( 'ena_sample_id', $ena_sample_id );
	if (@$runs) {
		return @$runs;
	}
	return undef;
}

sub store {
	my ( $self, $run, $update ) = @_;
	my $existing_record = $self->fetch_by_ena_sample_id( $run->ena_sample_id );

	if ( $existing_record && $update ) {
		$run->dbID( $existing_record->dbID );
		return $self->update($run);
	}

	if ( !$existing_record ) {
		$self->SUPER::store( $run, $update );
	}
}

1;
