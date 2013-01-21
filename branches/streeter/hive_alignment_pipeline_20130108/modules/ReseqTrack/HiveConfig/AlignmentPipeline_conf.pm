
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


package AlignmentPipeline_conf;

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

        'pipeline_name' => 'align',                     # name used by the beekeeper to prefix job names on the farm

        'reseqtrack_db'  => {
            -host => $self->o('host'),
            -port => 4197,
            -user => 'g1kro',
            -pass => '',
            #-dbname => 'nextgen_track', # set on the command line
        },

        'chunk_max_bases'    => 500000000,
        'type_fastq'    => 'FILTERED_FASTQ',
        'split_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/split/split',
        'bwa_exe' => '',
        'samtools_exe' => '',
        'squeeze_exe' => '',
        'gatk_dir' => '',
        'picard_dir' => '',

        'universal_branch_parameters_in' => {
          'branch_label' => 'label',
          'output_dir' => 'output_dir',
        },

    };
}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and populating it with Hive tables and procedures it also creates two pipeline-specific tables used by Runnables to communicate.

=cut

sub pipeline_create_commands {
    my ($self) = @_;
    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation

        $self->db_execute_command('pipeline_db', 'CREATE TABLE branch (child_branch_id int(10) unsigned NOT NULL AUTO_INCREMENT, parent_branch_id int(10) unsigned, child_branch_index int(10) unsigned, PRIMARY KEY (child_branch_id))'),
        $self->db_execute_command('pipeline_db', 'CREATE TABLE branch_meta_data (branch_meta_data_id int(10) unsigned NOT NULL AUTO_INCREMENT, branch_id int(10) unsigned NOT NULL, meta_key VARCHAR(50) NOT NULL, meta_value VARCHAR(1000) NOT NULL, is_active TINYINT(1), never_delete TINYINT(1), PRIMARY KEY (branch_meta_data_id))'),

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
        'universal_branch_parameters_in' => $self->o('universal_branch_parameters_in'),

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
    return [

        {   -logic_name    => 'get_runs',
            -module        => 'ReseqTrack::HiveProcess::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                type_branch => 'run',
                branch_parameters_out => {
                  split_fastq => 'FASTQ',
                },
            },
            -input_ids => [{'branch_id' => 0}],
            -flow_into => {
                '2->A' => [ 'split_fastq' ],   # will create a semaphored fan of jobs
                'A->1' => [ 'pipeline_done'  ],   # will create a semaphored funnel job to wait for the fan to complete
            },
        },

        {   -logic_name    => 'split_fastq',
            -module        => 'ReseqTrack::HiveProcess::SplitFastq',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                program_file => $self->o('split_exe'),
                type_fastq => $self->o('type_fastq'),
                max_bases => $self->o('chunk_max_bases'),
                branch_parameters_in => {
                    run_id => 'RUN_ID',
                },
                branch_parameters_out => {
                  split_fastq => 'FASTQ',
                },
            },
            -hive_capacity      => -1,  # turn off the reciprocal limiter
            -analysis_capacity  =>  4,  # use per-analysis limiter
            -flow_into => {
                '2->A' => [ 'bwa' ],   # will create a semaphored fan of jobs
                'A->1' => [ 'decide_merge_chunks'  ],   # will create a semaphored funnel job to wait for the fan to complete
            },
        },
        
        {   -logic_name => 'bwa',
            -module        => 'ReseqTrack::HiveProcess::BWA',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                program_file => $self->o('bwa_exe'),
                branch_parameters_in => {
                    fastq => {key => 'FASTQ', becomes_inactive => 1},
                    run_id => {key => 'RUN_ID', ascend => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -flow_into => {
                1 => [ 'sort_chunks'],
            },
        },

        {   -logic_name => 'sort_chunks',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'sort',
                branch_parameters_in => {
                    bam => {key => 'BAM', becomes_inactive => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
        },

        {   -logic_name => 'decide_merge_chunks',
            -module        => 'ReseqTrack::HiveProcess::FlowDecider',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                branch_parameters_in => {
                    file => {key => 'BAM', descend => 1},
                },
                flows_if_no_files => undef,
                flows_if_one_file => 1,
                flows_if_multiple_files => [2,1],
            },
            -flow_into => {
                '2->A' => [ 'merge_chunks'],
                'A->1' => [ 'mark_duplicates'],
            }
        },
        {   -logic_name => 'merge_chunks',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'merge',
                branch_parameters_in => {
                    bam => {key => 'BAM', descend => 1, becomes_inactive => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
        },
        {   -logic_name => 'mark_duplicates',
            -module        => 'ReseqTrack::HiveProcess::RunPicard',
            -meadow_type=> 'LOCAL',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                command => 'mark_duplicates',
                branch_parameters_in => {
                    bam => {key => 'BAM', becomes_inactive => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -flow_into => {
                1 => [ 'index_mrkdup_bam'],
            }
        },
        {   -logic_name => 'index_mrkdup_bam',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -meadow_type=> 'LOCAL',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'index',
                branch_parameters_in => {
                    bam => 'BAM',
                },
                branch_parameters_out => {
                  bai => 'BAI',
                },
            },
            -flow_into => {
                1 => [ 'realign_indels'],
            }
        },
        {   -logic_name => 'realign_indels',
            -module        => 'ReseqTrack::HiveProcess::RunBamImprovement',
            -meadow_type=> 'LOCAL',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'realign',
                reference => $self->o('reference'),
                gatk_dir => $self->o('gatk_dir'),
                jvm_args => '-Xmx2g',
                branch_parameters_in => {
                    bam => {key => 'BAM', becomes_inactive => 1},
                    bai => {key => 'BAI', becomes_inactive => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -flow_into => {
                1 => [ 'index_realigned_bam'],
            }
        },
        {   -logic_name => 'index_realigned_bam',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -meadow_type=> 'LOCAL',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'index',
                branch_parameters_in => {
                    bam => 'BAM',
                },
                branch_parameters_out => {
                  bai => 'BAI',
                },
            },
            -flow_into => {
                1 => [ 'recalibrate_bam'],
            }
        },
        {   -logic_name => 'recalibrate_bam',
            -module        => 'ReseqTrack::HiveProcess::RunBamImprovement',
            -meadow_type=> 'LOCAL',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'recalibrate',
                reference => $self->o('reference'),
                gatk_dir => $self->o('gatk_dir'),
                jvm_args => '-Xmx2g',
                branch_parameters_in => {
                    bam => {key => 'BAM', becomes_inactive => 1},
                    bai => {key => 'BAI', becomes_inactive => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -flow_into => {
                1 => [ 'calmd'],
            }
        },
        {   -logic_name => 'calmd',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'calmd',
                branch_parameters_in => {
                    bam => {key => 'BAM', becomes_inactive => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -flow_into => {
                1 => [ 'tag_strip'],
            }
        },
        {   -logic_name => 'tag_strip',
            -module        => 'ReseqTrack::HiveProcess::RunSqueezeBam',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                'program_file' => $self->o('squeeze_exe'),
                'rm_OQ_fields' => 1,
                'rm_tag_types' => ['XM:i', 'XG:i', 'XO:i'],
                'branch_parameters_in' => {
                    bam => {key => 'BAM', becomes_inactive => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
        },

        {   -logic_name => 'pipeline_done',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -module        => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;

