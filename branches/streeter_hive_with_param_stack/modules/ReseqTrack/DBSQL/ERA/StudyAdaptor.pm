package ReseqTrack::DBSQL::ERA::StudyAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::ERA::BaseEraAdaptor);
use ReseqTrack::Tools::Exception qw(throw);

sub table_name {
  return "study, cv_status, submission";
}

sub columns {
  return
"study.study_id, study.md5, study.study_type, study.study_title, study.study_alias, cv_status.status, submission.submission_id, to_char(submission.submission_date, 'YYYY-MM-DD HH24:MI') submission_date, study.ega_id";
}

sub where {
  return
"study.status_id = cv_status.status_id and study.submission_id = submission.submission_id";
}

sub object_from_hashref {
  my ( $self, $hashref ) = @_;
  throw("Can't create a ReseqTrack::Study from an undefined hashref")
    if ( !$hashref );

  my $study = ReseqTrack::Study->new(
    -source_id       => $hashref->{STUDY_ID},
    -md5             => $hashref->{MD5},
    -type            => $hashref->{STUDY_TYPE},
    -title           => $hashref->{STUDY_TITLE},
    -status          => $hashref->{STATUS},
    -submission_id   => $hashref->{SUBMISSION_ID},
    -submission_date => $hashref->{SUBMISSION_DATE},
    -study_alias     => $hashref->{STUDY_ALIAS}
  );
 $self->add_ega_id($study,$hashref);

  return $study;
}

sub internal_id_column {
  return "study_id";
}

sub fetch_by_study_id {
  my ( $self, $study_id ) = @_;
  return [ $self->fetch_by_dbID($study_id) ];
}

sub fetch_by_sample_id {
  my ( $self, $sample_id ) = @_;
  my $sql =
      "select distinct "
    . $self->columns
    . " from "
    . $self->table_name
    . ", run_sample, run, experiment";
  $sql .= " where " . $self->where if ( $self->where );
  $sql .= ' and run_sample.sample_id = ?';
  $sql .= ' and run_sample.run_id = run.run_id';
  $sql .= ' and run.experiment_id = experiment.experiment_id';
  $sql .= ' and experiment.study_id = study.study_id';

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

1;
