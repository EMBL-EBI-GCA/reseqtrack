package ReseqTrack::DBSQL::ERA::SampleAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::ERA::BaseEraAdaptor);
use ReseqTrack::Sample;
use ReseqTrack::Tools::Exception qw(throw);

sub table_name {
  return "sample, cv_status, submission";
}

sub columns {
  return
"sample.sample_id, cv_status.status, sample.md5, sample.center_name, sample.sample_alias, sample.tax_id, sample.scientific_name, sample.common_name, sample.anonymized_name, sample.individual_name, sample.sample_title, submission.submission_id, to_char(submission.submission_date, 'YYYY-MM-DD HH24:MI') submission_date, sample.ega_id, sample.biosample_id, sample.biosample_authority";
}

sub where {
  return
"sample.status_id = cv_status.status_id and sample.submission_id = submission.submission_id";
}

sub xml_column {
  return "sample.sample_xml";
}

sub object_from_hashref {
  my ( $self, $hashref ) = @_;
  throw("Can't create a ReseqTrack::Sample from an undefined hashref")
    if ( !$hashref );

  my $sample = ReseqTrack::Sample->new(
    -source_id           => $hashref->{SAMPLE_ID},
    -status              => $hashref->{STATUS},
    -md5                 => $hashref->{MD5},
    -center_name         => $hashref->{CENTER_NAME},
    -sample_alias        => $hashref->{SAMPLE_ALIAS},
    -tax_id              => $hashref->{TAX_ID},
    -scientific_name     => $hashref->{SCIENTIFIC_NAME},
    -common_name         => $hashref->{COMMON_NAME},
    -anonymized_name     => $hashref->{ANONYMIZED_NAME},
    -individual_name     => $hashref->{INDIVIDUAL_NAME},
    -sample_title        => $hashref->{SAMPLE_TITLE},
    -submission_id       => $hashref->{SUBMISSION_ID},
    -submission_date     => $hashref->{SUBMISSION_DATE},
    -biosample_id        => $hashref->{BIOSAMPLE_ID},
    -biosample_authority => $hashref->{BIOSAMPLE_AUTHORITY},
  );
  $self->add_ega_id( $sample, $hashref );

  return $sample;
}

sub attribute_tag {
  return 'SAMPLE_ATTRIBUTE';
}

sub internal_id_column {
  return "sample.sample_id";
}

sub fetch_by_study_id {
  my ( $self, $study_id ) = @_;
  my $sql = "select " . $self->columns . " from " . $self->table_name;
  $sql .= " where " . $self->where;
  $sql .= " and sample.sample_id in (";
  $sql .=
"  select EXTRACTVALUE(experiment_xml, '//EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR/\@accession') from";
  $sql .= "  experiment";
  $sql .= "  where experiment.study_id = ?";
  $sql .= ")";

  my @objects;
  my $sth = $self->prepare($sql);
  $sth->bind_param( 1, $study_id );
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
