=head1 NAME

 ReseqTrack::Hive::PipeConfig::Filerelease_conf

=head1 SYNOPSIS

  This is a pipeline for a file release pipeline.
  Files with a foreign host_id are moved from a dropbox to the project's file system.
  Messages from the pipeline get written to the attribute table of the ReseqTrack database.
  Rejected files get retried if they are updated in the ReseqTrack database or if their unix timestamp changes

  Pipeline must be seeded by the file table of a ReseqTrack database. (foreign files only)
  i.e. use the seeding module ReseqTrack::Hive::PipeSeed::ForeignFiles or something very similar to it

  Here is an example pipeline configuration to load using reseqtrack/scripts/pipeline/load_pipeline_from_conf.pl

[file release]
table_name=file
config_module=ReseqTrack::Hive::PipeConfig::FileRelease_conf
config_options=-file_move_module MyProjectModules::MoveFile

  Options that MUST be specified in the pipeline.config_options table/column of your ReseqTrack database:

      -file_move_module, A module derived from ReseqTrack::Hive::Process::FileRelease::Move
              This modules implements the derive_directory subrouine (i.e. a project-specific subroutine)

  Options that have defaults but you will often want to set them in your pipeline.cofig_options table/column:

      -seeding_module, (default is ReseqTrack::Hive::PipeSeed::ForeignFiles) override this with a project-specific module
      -seeding_options, hashref passed to the seeding module.  Default is {}.

      -checking_module, (default is ReseqTrack::Hive::Process::FileRelease::Checks)
      -hostname, (default is 1000genomes.ebi.ac.uk)

  Options that are required, but will be passed in by reseqtrack/scripts/init_pipeline.pl:

      -pipeline_db -host=???
      -pipeline_db -port=???
      -pipeline_db -user=???
      -dipeline_db -dbname=???
      -password
      -reseqtrack_db -host=???
      -reseqtrack_db -user=???
      -reseqtrack_db -port=???
      -reseqtrack_db -pass=???
      -reseqtrack_db -dbname=???

=cut


package ReseqTrack::Hive::PipeConfig::FileRelease_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'file_release',                     # name used by the beekeeper to prefix job names on the farm

        seeding_module => 'ReseqTrack::Hive::PipeSeed::ForeignFiles',
        seeding_options => {},

        checking_module  => 'ReseqTrack::Hive::Process::FileRelease::Checks',
        hostname => '1000genomes.ebi.ac.uk',

        derive_directory_options => {},

    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '200Mb' => { 'LSF' => '-C0 -M200 -q '.$self->o('lsf_queue').' -R"select[mem>200] rusage[mem=200]"' },
    };
}


sub pipeline_analyses {
    my ($self) = @_;

    my @analyses;
    push(@analyses, {
            -logic_name    => 'get_seeds',
            -module        => 'ReseqTrack::Hive::Process::SeedFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                seeding_module => $self->o('seeding_module'),
                seeding_options => $self->o('seeding_options'),
            },
            -flow_into => {
                2 => [ 'quick_checks' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'quick_checks',
            -module        => $self->o('checking_module'),
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
              check_class => 'quick',
              reseqtrack_options => {
                flows_non_factory => [1,2],
                flows_do_count => {
                  1 => '0',
                  2 => '1+',
                },
                flows_do_count_param => 'is_failed',
              },
            },
            -flow_into => {
                1 => [ 'slow_checks' ],
                2 => [ 'seed_complete' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'slow_checks',
            -module        => $self->o('checking_module'),
            -parameters    => {
              check_class => 'slow',
              reseqtrack_options => {
                flows_non_factory => [1,2],
                flows_do_count => {
                  1 => '0',
                  2 => '1+',
                },
                flows_do_count_param => 'is_failed',
              },
            },
            -flow_into => {
                1 => [ 'move_to_staging' ],
                2 => [ 'seed_complete' ],
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    push(@analyses, {
            -logic_name    => 'move_to_staging',
            -module        => $self->o('file_move_module'),
            -parameters    => {
                hostname => $self->o('hostname'),
                derive_directory_options => $self->o('derive_directory_options'),
            },
            -flow_into => {
                1 => ['seed_complete'],
            },
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
      });
    push(@analyses, {
            -logic_name    => 'seed_complete',
            -module        => 'ReseqTrack::Hive::Process::UpdateSeed',
            -parameters    => {
              delete_seeds  => '#expr(#is_failed# ? 0 : 1)expr#',
            },
            -meadow_type => 'LOCAL',
      });

    return \@analyses;
}

1;

