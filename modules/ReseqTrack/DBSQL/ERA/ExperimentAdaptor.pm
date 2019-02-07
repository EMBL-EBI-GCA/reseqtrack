package ReseqTrack::DBSQL::ERA::ExperimentAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::ERA::BaseEraAdaptor);
use ReseqTrack::Tools::Exception qw(throw);

use ReseqTrack::Experiment;

sub table_name {
  return "experiment, cv_status, submission";
}

sub columns {
  return
"experiment.experiment_id, experiment.study_id, cv_status.status, experiment.md5, experiment.center_name,
 experiment.experiment_alias, experiment.instrument_platform, experiment.instrument_model, experiment.library_layout,
 experiment.library_name, experiment.library_strategy, experiment.library_source, experiment.library_selection,
 experiment.paired_nominal_length, experiment.paired_nominal_sdev, submission.submission_id,
 to_char(submission.submission_date, 'YYYY-MM-DD HH24:MI') submission_date, experiment.ega_id,
 EXTRACTVALUE(experiment.experiment_xml, '//EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR/\@accession') sample_id,
 existsnode(experiment.experiment_xml,'//EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR/POOL') is_pool
";
}

sub where {
  return
"experiment.status_id = cv_status.status_id and experiment.submission_id = submission.submission_id";
}

sub xml_column {
  return "experiment.experiment_xml";
}

sub attribute_tag {
  return "EXPERIMENT_ATTRIBUTE";
}

sub object_from_hashref {
  my ( $self, $hashref ) = @_;
  throw("Can't create a ReseqTrack::Experiment from an undefined hashref")
    if ( !$hashref );

  my $exp = ReseqTrack::Experiment->new(
      -source_id             => $hashref->{EXPERIMENT_ID},
      -study_id              => $hashref->{STUDY_ID},
      -status                => $hashref->{STATUS},
      -md5                   => $hashref->{MD5},
      -center_name           => $hashref->{CENTER_NAME},
      -experiment_alias      => $hashref->{EXPERIMENT_ALIAS},
      -instrument_platform   => $hashref->{INSTRUMENT_PLATFORM},
      -instrument_model      => $hashref->{INSTRUMENT_MODEL},
      -library_layout        => $hashref->{LIBRARY_LAYOUT},
      -library_name          => $hashref->{LIBRARY_NAME},
      -library_strategy      => $hashref->{LIBRARY_STRATEGY},
      -library_source        => $hashref->{LIBRARY_SOURCE},
      -library_selection     => $hashref->{LIBRARY_SELECTION},
      -paired_nominal_length => $hashref->{PAIRED_NOMINAL_LENGTH},
      -paired_nominal_sdev   => $hashref->{PAIRED_NOMINAL_SDEV},
      -submission_id         => $hashref->{SUBMISSION_ID},
      -submission_date       => $hashref->{SUBMISSION_DATE},
      -sample_id             => $hashref->{SAMPLE_ID},
      -is_pool               => $hashref->{IS_POOL}
  );
  $self->add_ega_id( $exp, $hashref );

  return $exp;
}

sub internal_id_column {
  return "experiment.experiment_id";
}

sub fetch_by_study_id {
  my ( $self, $study_id ) = @_;
  return $self->fetch_by_column_name( "experiment.study_id", $study_id );
}

1;

