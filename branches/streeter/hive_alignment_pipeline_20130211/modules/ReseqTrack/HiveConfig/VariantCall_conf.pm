
=pod 

=head1 NAME

    Bio::EnsEMBL::Hive::PipeConfig::LongMult_conf;

=head1 SYNOPSIS

   # Example 1: specifying only the mandatory option (numbers to be multiplied are taken from defaults)
init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::LongMult_conf -password <mypass>

   # Example 2: specifying the mandatory options as well as overriding the default numbers to be multiplied:
init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::LongMult_conf -password <mypass> -first_mult 2344556 -second_mult 777666555

   # Example 3: do not re-create the database, just load another multiplicaton task into an existing one:
init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::LongMult_conf -job_topup -password <mypass> -first_mult 1111222233334444 -second_mult 38578377835


=head1 DESCRIPTION

    This is the PipeConfig file for the long multiplication pipeline example.
    The main point of this pipeline is to provide an example of how to write Hive Runnables and link them together into a pipeline.

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.


    The setting. Let's assume we are given two loooooong numbers to multiply. Reeeeally long.
    So long that they do not fit into registers of the CPU and should be multiplied digit-by-digit.
    For the purposes of this example we also assume this task is very computationally intensive and has to be done in parallel.

    The long multiplication pipeline consists of three "analyses" (types of tasks):  'start', 'part_multiply' and 'add_together'
    that we will be using to examplify various features of the Hive.

        * A 'start' job takes in two string parameters, 'a_multiplier' and 'b_multiplier',
          takes the second one apart into digits, finds what _different_ digits are there,
          creates several jobs of the 'part_multiply' analysis and one job of 'add_together' analysis.

        * A 'part_multiply' job takes in 'a_multiplier' and 'digit', multiplies them and records the result in 'intermediate_result' table.

        * An 'add_together' job waits for the first two analyses to complete,
          takes in 'a_multiplier', 'b_multiplier' and 'intermediate_result' table and produces the final result in 'final_result' table.

    Please see the implementation details in Runnable modules themselves.

=head1 CONTACT

    Please contact ehive-users@ebi.ac.uk mailing list with questions/suggestions.

=cut


package ReseqTrack::HiveConfig::VariantCall_conf;

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

        'pipeline_name' => 'vc',                     # name used by the beekeeper to prefix job names on the farm

        'reseqtrack_db'  => {
            -host => $self->o('host'),
            -port => 4197,
            -user => 'g1kro',
            -pass => '',
            #-dbname => 'nextgen_track', # set on the command line
        },

        'type_bam'    => 'BAM',
        'transpose_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/transpose_bam/transpose_bam',
        'samtools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/samtools',
        'bcftools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/bcftools/bcftools',
        'vcfutils_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/bcftools/vcfutils.pl',
        'bgzip_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
        'gatk_dir' => '/nfs/1000g-work/G1K/work/bin/gatk/dist/',
        'freebayes_exe' => '/nfs/1000g-work/G1K/work/bin/freebayes/bin/freebayes',

        call_by_gatk_options => {}, # use module defaults
        call_by_samtools_options => {}, # use module defaults
        call_by_freebayes_options => {}, # use module defaults

        'fai' => $self->o('reference') . '.fai',

        'final_label' => $self->o('pipeline_name'),
        'call_window_size' => 50000,


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

    my $sql_3 = "
    CREATE TABLE sequence_region (
      sequence_region_index int(10) unsigned NOT NULL,
      name VARCHAR(50) NOT NULL,
      length int(10) unsigned,
      PRIMARY KEY (sequence_region_index)
    )";

    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation

        $self->db_execute_command('pipeline_db', $sql_1),
        $self->db_execute_command('pipeline_db', $sql_2),
        $self->db_execute_command('pipeline_db', $sql_3),

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

        'root_output_dir' => $self->o('root_output_dir'),
        'reseqtrack_db' => $self->o('reseqtrack_db'),

        'universal_branch_parameters_in' => {
          'branch_label' => 'label',
          'branch_subdir' => {key => 'label', ascend => 1},
        },

    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '200Mb' => { 'LSF' => '-C0 -M200 -q production -R"select[mem>200] rusage[mem=200]"' },
            '500Mb' => { 'LSF' => '-C0 -M500 -q production -R"select[mem>500] rusage[mem=500]"' },
            '1Gb'   => { 'LSF' => '-C0 -M1000 -q production -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q production -R"select[mem>2000] rusage[mem=2000]"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q production -R"select[mem>4000] rusage[mem=4000]"' },
            '5Gb' => { 'LSF' => '-C0 -M5000 -q production -R"select[mem>5000] rusage[mem=5000]"' },
    };
}


=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.
                  Here it defines three analyses:

                    * 'start' with two jobs (multiply 'first_mult' by 'second_mult' and vice versa - to check the commutativity of multiplivation).
                      Each job will dataflow (create more jobs) via branch #2 into 'part_multiply' and via branch #1 into 'add_together'.

                    * 'part_multiply' initially without jobs (they will flow from 'start')

                    * 'add_together' initially without jobs (they will flow from 'start').
                       All 'add_together' jobs will wait for completion of 'part_multiply' jobs before their own execution (to ensure all data is available).

    There are two control modes in this pipeline:
        A. The default mode is to use the '2' and '1' dataflow rules from 'start' analysis and a -wait_for rule in 'add_together' analysis for analysis-wide synchronization.
        B. The semaphored mode is to use '2->A' and 'A->1' semaphored dataflow rules from 'start' instead, and comment out the analysis-wide -wait_for rule, relying on semaphores.

=cut

sub pipeline_analyses {
    my ($self) = @_;

    my $call_by_samtools = 1;
    my $call_by_gatk = 1;
    my $call_by_freebayes = 1;

    my @call_processes;
    push(@call_processes, 'call_by_samtools') if $call_by_samtools;
    push(@call_processes, 'call_by_gatk') if $call_by_gatk;
    push(@call_processes, 'call_by_freebayes') if $call_by_freebayes;

    my @merge_processes;
    push(@merge_processes, 'merge_samtools') if $call_by_samtools;
    push(@merge_processes, 'merge_gatk') if $call_by_gatk;
    push(@merge_processes, 'merge_freebayes') if $call_by_freebayes;

    my @analyses;
    push(@analyses, {
            -logic_name    => 'import_genome',
            -module        => 'ReseqTrack::HiveProcess::RegionsFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                command => 'load_db',
                fai => $self->o('fai'),
            },
            -input_ids => [{}],
      });
    push(@analyses, {
            -logic_name    => 'callgroups_factory',
            -module        => 'ReseqTrack::HiveProcess::JobFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                branch_parameters_out => {
                  data_value => 'callgroup',
                }
            },
            -input_ids => [{data_values => [split(',', $self->o('callgroup'))]}],
            -flow_into => {
                2 => [ 'find_source_bams' ],   # will create a semaphored fan of jobs
            },
      });
    push(@analyses, {
            -logic_name    => 'find_source_bams',
            -module        => 'ReseqTrack::HiveProcess::ImportCollection',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                collection_type => $self->o('type_bam'),
                branch_parameters_in => {
                    collection_name => 'callgroup',
                },
                branch_parameters_out => {
                  name => 'SOURCE_BAM',
                },
            },
            -flow_into => {
                1 => [ 'regions_factory_1' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'regions_factory_1',
            -module        => 'ReseqTrack::HiveProcess::RegionsFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                command => 'seq_region_factory',
                branch_parameters_out => {
                  seq_index_start => 'seq_index_start',
                  seq_index_end => 'seq_index_end',
                },
                min_bases => 50000000,
            },
            -wait_for => ['import_genome'],
            -flow_into => {
                '2->A' => [ 'transpose_bam' ],
                'A->1' => [ 'dummy_process' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'transpose_bam',
            -module        => 'ReseqTrack::HiveProcess::TransposeBam',
            -parameters    => {
                program_file => $self->o('transpose_exe'),
                create_index => 1,
                branch_parameters_in => {
                    bams => {key => 'SOURCE_BAM', ascend => 1},
                    seq_index_start => 'seq_index_start',
                    seq_index_end => 'seq_index_end',
                },
                branch_parameters_out => {
                  bam => 'BAM',
                  bai => 'BAI',
                },
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  4,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                '1->A' => [ 'regions_factory_2' ],
                'A->1' => [ 'delete_bam' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'regions_factory_2',
            -module        => 'ReseqTrack::HiveProcess::RegionsFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                command => 'window_factory',
                branch_parameters_in => {
                    seq_index_start => 'seq_index_start',
                    seq_index_end => 'seq_index_end',
                },
                branch_parameters_out => {
                  seq_index => 'seq_index',
                  region_start => 'region_start',
                  region_end => 'region_end',
                },
                max_bases => $self->o('call_window_size'),
            },
            -flow_into => {
                2 => [ @call_processes ],
            },
      });
    push(@analyses, {
            -logic_name    => 'delete_bam',
            -module        => 'ReseqTrack::HiveProcess::BranchableProcess',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                branch_parameters_in => {
                    bams => {key => 'BAM', inactivate => 1},
                    bais => {key => 'BAI', inactivate => 1},
                },
            },
      });
    push(@analyses, {
            -logic_name    => 'dummy_process',
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -flow_into => {
                1 => [ @merge_processes ],
            }
      });
    if ($call_by_samtools) {
      push(@analyses, {
            -logic_name    => 'call_by_samtools',
            -module        => 'ReseqTrack::HiveProcess::RunVariantCall::CallBySamtools',
            -parameters    => {
                reference => $self->o('reference'),
                samtools => $self->o('samtools_exe'),
                bcftools => $self->o('bcftools_exe'),
                vcfutils => $self->o('vcfutils_exe'),
                bgzip => $self->o('bgzip_exe'),
                options => $self->o('call_by_samtools_options'),
                branch_parameters_in => {
                    seq_index => 'seq_index',
                    region_start => 'region_start',
                    region_end => 'region_end',
                    bam => {key => 'BAM', ascend => 1},
                },
                branch_parameters_out => {
                  vcf => 'samtools_vcf'
                },
            },
            -rc_name => '500Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
        });
      push(@analyses, {
            -logic_name    => 'merge_samtools',
            -module        => 'ReseqTrack::HiveProcess::MergeVcf',
            -parameters    => {
                bgzip => $self->o('bgzip_exe'),
                branch_parameters_in => {
                  vcf => {key => 'samtools_vcf', inactivate => 1, descend => 1},
                },
                branch_parameters_out => {
                  vcf => 'samtools_vcf'
                },
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
        });
    }
    if ($call_by_gatk) {
      push(@analyses, {
            -logic_name    => 'call_by_gatk',
            -module        => 'ReseqTrack::HiveProcess::RunVariantCall::CallByGATK',
            -parameters    => {
                reference => $self->o('reference'),
                gatk_dir => $self->o('gatk_dir'),
                options => $self->o('call_by_gatk_options'),
                branch_parameters_in => {
                    seq_index => 'seq_index',
                    region_start => 'region_start',
                    region_end => 'region_end',
                    bam => {key => 'BAM', ascend => 1},
                },
                branch_parameters_out => {
                  vcf => 'gatk_vcf'
                },
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
        });
      push(@analyses, {
            -logic_name    => 'merge_gatk',
            -module        => 'ReseqTrack::HiveProcess::MergeVcf',
            -parameters    => {
                bgzip => $self->o('bgzip_exe'),
                branch_parameters_in => {
                  vcf => {key => 'gatk_vcf', inactivate => 1, descend => 1},
                },
                branch_parameters_out => {
                  vcf => 'gatk_vcf'
                },
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
        });
    }
    if ($call_by_freebayes) {
      push(@analyses, {
            -logic_name    => 'call_by_freebayes',
            -module        => 'ReseqTrack::HiveProcess::RunVariantCall::CallByFreebayes',
            -parameters    => {
                reference => $self->o('reference'),
                freebayes => $self->o('freebayes_exe'),
                bgzip => $self->o('bgzip_exe'),
                options => $self->o('call_by_freebayes_options'),
                branch_parameters_in => {
                    seq_index => 'seq_index',
                    region_start => 'region_start',
                    region_end => 'region_end',
                    bam => {key => 'BAM', ascend => 1},
                },
                branch_parameters_out => {
                  vcf => 'freebayes_vcf'
                },
            },
            -rc_name => '500Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
        });
      push(@analyses, {
            -logic_name    => 'merge_freebayes',
            -module        => 'ReseqTrack::HiveProcess::MergeVcf',
            -parameters    => {
                bgzip => $self->o('bgzip_exe'),
                branch_parameters_in => {
                  vcf => {key => 'freebayes_vcf', inactivate => 1, descend => 1},
                },
                branch_parameters_out => {
                  vcf => 'freebayes_vcf'
                },
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
        });
    }

    return \@analyses;
}

1;

