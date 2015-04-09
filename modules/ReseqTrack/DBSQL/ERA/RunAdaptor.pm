package ReseqTrack::DBSQL::ERA::RunAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::ERA::BaseEraAdaptor);
use ReseqTrack::Run;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::AttributeUtils qw(create_attribute_for_object);

sub table_name {
    return "run, experiment,  cv_status, submission";
}

sub columns {
    return
"run.run_id, run.experiment_id, EXTRACTVALUE(experiment_xml, '//EXPERIMENT/DESIGN/SAMPLE_DESCRIPTOR/\@accession') sample_id, run.run_alias, cv_status.status, run.md5, run.center_name, run.run_center_name, run.instrument_platform, run.instrument_model, submission.submission_id, to_char(submission.submission_date, 'YYYY-MM-DD HH24:MI')   submission_date,run.ega_id,replace(EXTRACTVALUE(run.run_xml, '/RUN_SET/RUN/\@run_date'),'T',' ') run_date";
}

sub where {
    return
"run.status_id = cv_status.status_id and run.submission_id = submission.submission_id and experiment.experiment_id = run.experiment_id";
}

sub object_from_hashref {
    my ( $self, $hashref ) = @_;
    throw("Can't create a ReseqTrack::Run from an undefined hashref")
      if ( !$hashref );

    my $run = ReseqTrack::Run->new(
        -source_id           => $hashref->{RUN_ID},
        -md5                 => $hashref->{MD5},
        -experiment_id       => $hashref->{EXPERIMENT_ID},
        -sample_id           => $hashref->{SAMPLE_ID},
        -run_alias           => $hashref->{RUN_ALIAS},
        -status              => $hashref->{STATUS},
        -center_name         => $hashref->{CENTER_NAME},
        -run_center_name     => $hashref->{RUN_CENTER_NAME},
        -instrument_platform => $hashref->{INSTRUMENT_PLATFORM},
        -instrument_model    => $hashref->{INSTRUMENT_MODEL},
        -submission_id       => $hashref->{SUBMISSION_ID},
        -submission_date     => $hashref->{SUBMISSION_DATE},
        -run_date            => $hashref->{RUN_DATE},
    );
    $self->add_ega_id( $run, $hashref );
    $self->add_attribute( $run, $hashref, 'RUN_DATE', 'RUN_DATE' ) if ($hashref->{RUN_DATE});

    return $run;
}

sub add_ega_id {
    my ( $self, $object, $hashref ) = @_;

    if ( $hashref->{RUN_DATE} ) {
        my $attr =
          create_attribute_for_object( $object, 'RUN_DATE',
            $hashref->{RUN_DATE} );
        $object->attributes( [$attr] );
    }
}

sub get_run_stats {
    my ( $self, $run_id ) = @_;

    my $sql =
"select spot_count as read_count, base_count from run_stats where run_id = ?";
    my $sth = $self->prepare($sql);
    $sth->bind_param( 1, $run_id );
    $sth->execute;
    if ($@) {
        throw("Problem running $sql $@");
    }
    my $hashref = $sth->fetchrow_hashref;
    $sth->finish();
    return $hashref;
}

sub get_run_process {
    my ( $self, $run_id ) = @_;

    my $sql =
      "select fastq_date, fastq_error from run_process where run_id = ?";
    my $sth = $self->prepare($sql);
    $sth->bind_param( 1, $run_id );
    $sth->execute;
    if ($@) {
        throw("Problem running $sql $@");
    }
    my $hashref = $sth->fetchrow_hashref;
    $sth->finish();
    return $hashref;
}

sub internal_id_column {
    return "run.run_id";
}

sub xml_column {
    return "run_xml";
}

sub attribute_tag {
    return "RUN_ATTRIBUTE";
}

sub fetch_by_study_id {
    my ( $self, $study_id ) = @_;
    my $sql = "select " . $self->columns . " from " . $self->table_name;
    $sql .= " where " . $self->where;
    $sql .= " and experiment.study_id =  ?";

    print STDERR $sql.$/;

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
