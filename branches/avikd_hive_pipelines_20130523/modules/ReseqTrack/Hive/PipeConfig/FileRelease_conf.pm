
package ReseqTrack::Hive::PipeConfig::FileRelease_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');


=head2 default_options

    Description : Implements default_options() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.
                  In addition to the standard things it defines two options, 'first_mult' and 'second_mult' that are supposed to contain the long numbers to be multiplied.

=cut

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'file_release',                     # name used by the beekeeper to prefix job names on the farm

        checking_module  => 'ReseqTrack::Hive::Process::FileRelease::Checks',
        hostname => '1000genomes.ebi.ac.uk',


        use_label_management => 0,
        use_reseqtrack_file_table => 1,

    };
}

sub forbid_duplicate_file_id {return 1;}


sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '200Mb' => { 'LSF' => '-C0 -M200 -q production -R"select[mem>200] rusage[mem=200]"' },
    };
}


sub pipeline_analyses {
    my ($self) = @_;

    my @analyses;
    push(@analyses, {
            -logic_name    => 'foreign_files_factory',
            -module        => 'ReseqTrack::Hive::Process::FileRelease::ForeignFilesFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -input_ids => [{}],
            -flow_into => {
                '2->A' => [ 'quick_checks' ],   # will create a semaphored fan of jobs
                'A->1' => [ 'remove_lock' ],   # will create a semaphored fan of jobs
            },
      });
    push(@analyses, {
            -logic_name    => 'quick_checks',
            -module        => $self->o('checking_module'),
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
              check_class => 'quick',
            },
            -flow_into => {
                1 => [ 'slow_checks' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'slow_checks',
            -module        => $self->o('checking_module'),
            -parameters    => {
              check_class => 'slow',
            },
            -flow_into => {
                1 => [ 'move_to_staging' ],
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
            },
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
      });
    push(@analyses, {
            -logic_name    => 'remove_lock',
            -module        => 'ReseqTrack::Hive::Process::FileRelease::RemoveLock',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
      });

    return \@analyses;
}

1;

