
=head1 NAME

ReseqTrack::Tools::UpdateMetaData

=head1 SYNOPSIS

Object to update and sanity check the run meta info tables

=head1 Example


=cut

package ReseqTrack::Tools::UpdateMetaData;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::AttributeUtils
  qw(remove_outdated_attributes create_attribute_for_object);
use ReseqTrack::Tools::Metadata::EnaReadInfo;
use feature 'switch';
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception;

sub new {
    my ( $class, @args ) = @_;
    my $self = {};
    bless $self, $class;
    my (
        $era_db,   $dcc_db,          $add_ins,      $verbose,
        $load_new, $update_existing, $target_types, $target_studies,
        $use_rsm,  $run_stat_add_in, $log_fh
      )
      = rearrange(
        [
            qw(ERA_DB DCC_DB ADD_INS VERBOSE LOAD_NEW UPDATE_EXISTING TARGET_TYPES TARGET_STUDIES USE_DEFAULT_RSM RUN_STAT_ADD_IN LOG_FH )
        ],
        @args
      );

    $self->era_db($era_db);
    $self->dcc_db($dcc_db);
    $self->add_ins($add_ins);
    $self->verbose($verbose);
    $self->target_types($target_types);
    $self->target_studies($target_studies);
    $self->use_default_rsm($use_rsm);
    $self->run_stat_add_in($run_stat_add_in);
    $self->log_fh($log_fh);

    return $self;
}

sub load_new_studies {
    my ( $self, @studies_to_add ) = @_;
    for my $study_id (@studies_to_add) {
        $self->log("Adding study $study_id");
        $self->load_type_by_study_id( 'study', $study_id, 0, 1, 0 );
    }
}

sub update_from_era {
    my ( $self, $load_new, $update_existing, $force_update ) = @_;

    my $studies        = $self->dcc_db->get_StudyAdaptor()->fetch_all();
    my $target_studies = $self->target_studies;
    my $types          = $self->target_types;

  STUDY: for my $study (@$studies) {
        my @study_match = ( grep { $_ eq $study->source_id } @$target_studies );
        if ( @$target_studies && !@study_match ) {
            $self->log( "Skipping " . $study->source_id );
            next STUDY;
        }

        $self->log( "Checking " . $study->source_id );
      TYPE: for my $type (@$types) {
            $self->log("...for $type") if $self->verbose();
            $self->load_type_by_study_id( $type, $study->source_id(),
                $update_existing, $load_new, $force_update );
        }
    }
}

sub load_type_by_study_id {
    my ( $self, $type, $study_id, $update_existing, $load_new, $force_update )
      = @_;

    my $era_db   = $self->era_db;
    my $reseq_db = $self->dcc_db;
    my $verbose  = $self->verbose;
    my $add_ins  = $self->add_ins;

    my $era_adaptor   = adaptors( $era_db,   $type );
    my $reseq_adaptor = adaptors( $reseq_db, $type );
    my $fk_handler_sub        = fk_handlers($type);
    my $on_update_handler_sub = on_update_handlers($type);
    my $allow_attributes      = have_attributes($type);
    my $run_stat_add_in       = $self->run_stat_add_in;

    my $objects = $era_adaptor->fetch_by_study_id($study_id);
    die("Undef returned when fetching by study_id - check study ID is valid")
      if ( !defined $objects );

    my $stored_count  = 0;
    my $checked_count = scalar(@$objects);

    for my $object (@$objects) {
        my $current_record =
          $reseq_adaptor->fetch_by_source_id( $object->source_id );

        $fk_handler_sub->( $object, $reseq_db );
        my $run_stat_update;
        if ($run_stat_add_in) {
            $run_stat_update =
              $run_stat_add_in->check( $object, $current_record );
        }

        my $do_update = (
            $force_update || ( $update_existing
                && defined $current_record
                && changed( $object, $current_record ) )
        );

        my $do_load = ( $load_new && !defined $current_record );

        if ( $do_load || $do_update || $run_stat_update ) {
            $era_adaptor->attach_attributes($object) if ($allow_attributes);
            for (@$add_ins) {
                $_->check( $object, $current_record );
            }
            $self->log( "Storing $type " . $object->source_id() );

            $object->dbID( $current_record->dbID() ) if ($current_record);
            $reseq_adaptor->store( $object, 1 );

            if ($allow_attributes) {
                $reseq_adaptor->store_attributes( $object, 1 );
                remove_outdated_attributes( $object, $current_record,
                    $reseq_db->get_AttributeAdaptor )
                  if ($current_record);
            }
            $stored_count++;
        }
        if ($do_update) {
            $on_update_handler_sub->( $object, $reseq_db );
        }
    }

    $self->log("$type checked $checked_count stored $stored_count");
}

sub dcc_db {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{dcc_db} = $arg;
    }
    return $self->{dcc_db};
}

sub era_db {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{era_db} = $arg;
    }
    return $self->{era_db};
}

sub log_fh {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{log_fh} = $arg;
    }
    return $self->{log_fh};
}

sub add_ins {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{add_ins} = $arg;
    }
    return $self->{add_ins};
}

sub use_default_rsm {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{use_rsm} = $arg;
    }
    return $self->{use_rsm};
}

sub run_stat_add_in {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{rsm} = $arg;
    }
    if ( !$self->{rsm} && $self->use_default_rsm ) {
        $self->{rsm} = ReseqTrack::Tools::Metadata::EnaReadInfo->new(
            -era_db   => $self->era_db,
            -reseq_db => $self->dcc_db
        );
    }

    return $self->{rsm};
}

sub target_types {
    my ( $self, $target_types ) = @_;

# the order of types matters - e.g. experiments reference studies, runs reference samples and experiments
    my @types = qw(study sample experiment run);

    if ( $target_types && @$target_types ) {
        my %valid_targets;
        map { $valid_targets{$_} = 1 } @$target_types;
        @types = grep { $valid_targets{$_} } @types;
        $self->{target_types} = \@types;
    }

    if ( !$self->{target_types} ) {
        $self->{target_types} = \@types;
    }

    return $self->{target_types};
}

sub target_studies {
    my ( $self, $target_studies ) = @_;
    if ($target_studies) {
        my @ts = sort { $a cmp $b } @$target_studies;
        $self->{target_studies} = \@ts;
    }
    return $self->{target_studies};
}

sub changed {
    my ( $archive, $copy ) = @_;
    return ( $archive->md5 ne $copy->md5 || $archive->status ne $copy->status );
}

sub fk_handlers {
    my ($t) = @_;
    given ($t) {
        when ('study') {
            return sub { }
        }
        when ('sample') {
            return sub { }
        }
        when ('experiment') {
            return sub {
                my ( $e, $db ) = @_;

                throw( "No study ID set for experiment " . $e->source_id() )
                  unless ( $e->study_id );

                my $st_a = $db->get_StudyAdaptor();
                my $stid = $e->study_id();
                my $st   = $st_a->fetch_by_source_id($stid);
                throw( "Could not find study $stid for experiment "
                      . $e->source_id() )
                  unless $st;
                $e->study_id( $st->dbID() );

              }
        }
        when ('run') {
            return sub {
                my ( $r, $db ) = @_;

                throw( "No experiment ID set for run " . $r->source_id() )
                  unless ( $r->experiment_id );
                throw( "No sample ID set for run " . $r->source_id() )
                  unless ( $r->sample_id );

                my $ea  = $db->get_ExperimentAdaptor();
                my $eid = $r->experiment_id();
                my $e   = $ea->fetch_by_source_id($eid);
                throw( "Could not find experiment $eid for run "
                      . $r->source_id() )
                  unless $e;
                $r->experiment_id( $e->dbID() );

                my $sa_a = $db->get_SampleAdaptor();
                my $said = $r->sample_id();
                my $sa   = $sa_a->fetch_by_source_id($said);
                throw(
                    "Could not find sample $said for run " . $r->source_id() )
                  unless $sa;
                $r->sample_id( $sa->dbID() );
              }
        }
    }
}

# the db adaptors required for each data type
sub adaptors {
    my ( $db, $t ) = @_;
    given ($t) {
        when ('study')      { return $db->get_StudyAdaptor() }
        when ('sample')     { return $db->get_SampleAdaptor() }
        when ('experiment') { return $db->get_ExperimentAdaptor() }
        when ('run')        { return $db->get_RunAdaptor() }
    }
}

# which data types can use the attribute table each data type
sub have_attributes {
    my ($t) = @_;
    given ($t) {
        when ('study')      { return 0 }
        when ('sample')     { return 1 }
        when ('experiment') { return 1 }
        when ('run')        { return 1 }
    }
}

sub on_update_handlers {
    my ($t) = @_;
    given ($t) {
        when ('study') {
            return sub { }
        }
        when ('sample') {
            return sub { }
        }
        when ('experiment') {
            return sub {
                my ( $e, $db ) = @_;
                my $ra = $db->get_RunAdaptor();
                my $runs = $ra->fetch_by_experiment_id($e->dbID());
                for my $r (@$runs){
                  $r->md5('force_refresh');
                  $ra->update($r);
                }                
              }
        }
        when ('run') {
            return sub { }
        }
    }
}

sub log {
    my ( $self, $message ) = @_;
    my $log_fh = $self->log_fh;

    if ($log_fh) {
        print $log_fh $message . $/;
    }
}

sub report {
    my ($self) = @_;
    map { $_->report() } @{ $self->add_ins };
}

1;
