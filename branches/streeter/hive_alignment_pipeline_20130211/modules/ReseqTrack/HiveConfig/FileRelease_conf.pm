
package ReseqTrack::HiveConfig::FileRelease_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');  # All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.
                  In addition to the standard things it defines two options, 'first_mult' and 'second_mult' that are supposed to contain the long numbers to be multiplied.

=cut

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'file_release',                     # name used by the beekeeper to prefix job names on the farm

        'reseqtrack_db'  => {
            -host => $self->o('host'),
            -port => 4197,
            -user => 'g1krw',
            #-pass => '', # set on the command line
            #-dbname => 'file_release', # set on the command line
        },

        checking_module  => 'ReseqTrack::HiveProcess::FileReleaseChecks',

    };
}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and populating it with Hive tables and procedures it also creates two pipeline-specific tables used by Runnables to communicate.

=cut

sub pipeline_create_commands {
    my ($self) = @_;

    my $sql_1 = '
    CREATE TABLE branch (
      branch_id int(10) unsigned NOT NULL AUTO_INCREMENT,
      parent_branch_id int(10) unsigned,
      sibling_index int(10) unsigned,
      PRIMARY KEY (branch_id)
    )';

    my $sql_2 = "
    CREATE TABLE branch_data (
      branch_data_id int(10) unsigned NOT NULL AUTO_INCREMENT,
      branch_id int(10) unsigned NOT NULL,
      data_key VARCHAR(50) NOT NULL,
      data_value VARCHAR(1000) NOT NULL,
      is_active TINYINT(1),
      PRIMARY KEY (branch_data_id)
    )";

    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation

        $self->db_execute_command('pipeline_db', $sql_1),
        $self->db_execute_command('pipeline_db', $sql_2),

    ];
}


=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.

=cut

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class

        'reseqtrack_db' => $self->o('reseqtrack_db'),

        'universal_branch_parameters_in' => {
          'dropbox_filename' => 'filename',
          'db_md5' => 'md5',
          'db_size' => 'size',
          'db_file_id' => 'file_id',
        },

    };
}

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
            -module        => 'ReseqTrack::HiveProcess::ForeignFilesFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -input_ids => [{}],
            -flow_into => {
                '2->A' => [ 'quick_checks' ],   # will create a semaphored fan of jobs
                'A->1' => [ 'foreign_files_factory'  ],   # will create a semaphored funnel job to wait for the fan to complete
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
            -module        => 'ReseqTrack::HiveProcess::FileReleaseMove',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                branch_parameters_in => {
                    branch_timestamp => (key => 'factory_timestamp', ascend => 1),
                },
            },
      });

    return \@analyses;
}

1;

