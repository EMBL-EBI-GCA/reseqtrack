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
    experiment_source_id             => sub { $e->experiment_source_id(@_) },
    submission_id         => sub { $e->submission_id(@_) },
    submission_date       => sub { $e->submission_date(@_) },
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
  return pop @{ $self->fetch_by_column_name( "experiment_source_id", $source_id ) };
}

sub fetch_by_sample_id {
  my ( $self, $sample_id ) = @_;

  my $sql = join( ' ',
    'select distinct',
    $self->columns,
    'from',
    $self->table_name,
    ', run','where',
    'run.sample_id = ?',
    'and experiment.experiment_id = run.experiment_id' );

  my @objects;
  my $sth = $self->prepare($sql);
  $sth->bind_param( 1, $sample_id );

  eval { $sth->execute; };
  if ($@) {
    throw("Problem running $sql $@");
  }
  while ( my $rowHashref = $sth->fetchrow_hashref ) {
    my $object = $self->object_from_hashref($rowHashref) if ($rowHashref);
    push( @objects, $object );
  }
  $sth->finish;
  return \@objects;
}

sub fetch_by_sample_id_and_study_id {
  my ( $self, $sample_id, $study_id) = @_;

  my $sql = join( ' ',
    'select distinct',
    $self->columns,
    'from',
    $self->table_name,
    ', run','where',
    'and run.sample_id = ?',
    'and experiment.experiment_id = run.experiment_id',
    'and experiment.study_id = ?' 
    );

  my @objects;
  my $sth = $self->prepare($sql);
  $sth->bind_param( 1, $sample_id );
  $sth->bind_param( 2, $study_id );
  
  eval { $sth->execute; };
  if ($@) {
    throw("Problem running $sql $@");
  }
  while ( my $rowHashref = $sth->fetchrow_hashref ) {
    my $object = $self->object_from_hashref($rowHashref) if ($rowHashref);
    push( @objects, $object );
  }
  $sth->finish;
  return \@objects;
}

sub store {
  my ( $self, $experiment, $update ) = @_;
  my $existing_record = $self->fetch_by_dbID( $experiment->dbID )
    if ( $experiment->dbID );

  if ( $existing_record && $update ) {
    $experiment->dbID( $existing_record->dbID );
    return $self->update($experiment);
  }

  if ( !$existing_record ) {
    $self->SUPER::store( $experiment, $update );
  }
}

1;
