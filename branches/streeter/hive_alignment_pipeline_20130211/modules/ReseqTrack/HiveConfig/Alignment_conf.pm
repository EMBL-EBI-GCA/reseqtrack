
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


package ReseqTrack::HiveConfig::Alignment_conf;

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
        'validate_bam_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/validate_bam/validate_bam',
        'bwa_exe' => '/nfs/1000g-work/G1K/work/bin/bwa/bwa',
        'samtools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/samtools',
        'squeeze_exe' => '/nfs/1000g-work/G1K/work/bin/bamUtil/bin/bam',
        'gatk_dir' => '/nfs/1000g-work/G1K/work/bin/gatk/dist/',
        'picard_dir' => '/nfs/1000g-work/G1K/work/bin/picard',
        'known_indels_vcf' => '',
        'known_snps_vcf' => '',

        'reference_uri' => $self->o('reference'),

        'final_label' => $self->o('pipeline_name'),


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

    my $recalibrate_level = 0;
    my $realign_knowns_only = 0;

    die("recalibrate level must not be >2") if $recalibrate_level >2;

    # These are defaults for realign_knowns_only==0 && recalibrate_level==2
    my $process_run_level = 0; #
    my $do_index_for_recalibrate = 1; #
    my $flow_split_fastq = { 2 => ['bwa']}; #
    my $flow_decide_merge_chunks = undef; #
    my $flow_decide_merge_libraries = {'2->A' => ['merge_libraries'], 'A->1' => ['realign']}; #
    my $flow_realign = { 1 => ['index_for_recalibrate']}; #
    my $flow_recalibrate = { 1 => ['calmd']}; #
    my $flow_tag_strip = { 1 => ['reheader']}; #

    if (! $realign_knowns_only && $recalibrate_level == 0) {
      $do_index_for_recalibrate = 0;
      $flow_realign = { 1 => ['calmd']};
    }
    elsif (! $realign_knowns_only && $recalibrate_level == 1) {
      $process_run_level = 1;
      $do_index_for_recalibrate = 0;
      $flow_split_fastq = {'2->A' => ['bwa'], 'A->1' => ['decide_merge_chunks']};
      $flow_decide_merge_chunks = {'2->A' => ['merge_chunks'], 'A->1' => ['recalibrate']};
      $flow_realign = { 1 => ['calmd']};
      $flow_recalibrate = {};
    }
    elsif ($realign_knowns_only && $recalibrate_level == 1) {
      $process_run_level = 1;
      $do_index_for_recalibrate = 1;
      $flow_split_fastq = {'2->A' => [ 'bwa' ], 'A->1' => ['decide_merge_chunks'], };
      $flow_decide_merge_chunks = {'2->A' => ['merge_chunks'], 'A->1' => ['realign']};
      $flow_decide_merge_libraries = {2 => ['merge_libraries']};
      $flow_decide_merge_libraries = {'2->A' => ['merge_libraries'], 'A->1' => ['reheader']};
      $flow_tag_strip = {};
    }
    elsif ($realign_knowns_only && $recalibrate_level == 0) {
      $process_run_level = 1;
      $do_index_for_recalibrate = 0;
      $flow_split_fastq = {'2->A' => [ 'bwa' ], 'A->1' => ['decide_merge_chunks'], };
      $flow_decide_merge_chunks = {'2->A' => ['merge_chunks'], 'A->1' => ['realign']};
      $flow_decide_merge_libraries = {'2->A' => ['merge_libraries'], 'A->1' => ['reheader']};
      $flow_realign = { 1 => ['calmd']};
      $flow_tag_strip = {};
    }
    elsif ($realign_knowns_only && $recalibrate_level == 2) {
      die("pipeline not supported: realign_knowns_only==1 && recalibrate_level==2");
    }



    my @analyses;
    push(@analyses, {
            -logic_name    => 'studies_factory',
            -module        => 'ReseqTrack::HiveProcess::JobFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                branching_system => 1,
                branch_parameters_out => {
                  data_value => 'STUDY_ID',
                }
            },
            -input_ids => [{data_values => [split(',', $self->o('study_id'))]}],
            -flow_into => {
                2 => [ 'samples_factory' ],   # will create a semaphored fan of jobs
            },
      });
    push(@analyses, {
            -logic_name    => 'samples_factory',
            -module        => 'ReseqTrack::HiveProcess::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                type_branch => 'sample',
                branch_parameters_in => {
                    branch_study_id => 'STUDY_ID',
                }
            },
            -flow_into => {
                2 => [ 'libraries_factory' ],   # will create a semaphored fan of jobs
            },
      });
    push(@analyses, {
            -logic_name    => 'libraries_factory',
            -module        => 'ReseqTrack::HiveProcess::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                type_branch => 'library',
                branch_parameters_in => {
                    branch_sample_id => 'SAMPLE_ID',
                },
            },
            -flow_into => {
                '2->A' => [ 'runs_factory' ],   # will create a semaphored fan of jobs
                'A->1' => [ 'decide_merge_libraries'  ],   # will create a semaphored funnel job to wait for the fan to complete
            },
      });
    push(@analyses, {
            -logic_name    => 'runs_factory',
            -module        => 'ReseqTrack::HiveProcess::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                type_branch => 'run',
                branch_parameters_in => {
                    branch_sample_id => {key => 'SAMPLE_ID', ascend => 1},
                    branch_library_name => 'LIBRARY_NAME',
                },
            },
            -flow_into => {
                '2->A' => [ 'split_fastq' ],   # will create a semaphored fan of jobs
                'A->1' => [ 'decide_merge_runs'  ],   # will create a semaphored funnel job to wait for the fan to complete
            },
      });
    push(@analyses, {
            -logic_name    => 'split_fastq',
            -module        => 'ReseqTrack::HiveProcess::SplitFastq',
            -parameters    => {
                program_file => $self->o('split_exe'),
                type_fastq => $self->o('type_fastq'),
                max_bases => $self->o('chunk_max_bases'),
                branch_parameters_in => {
                    run_id => 'RUN_ID',
                },
                branch_parameters_out => {
                  split_fastq => 'FASTQ',
                  source_fastq => 'SOURCE_FASTQ',
                },
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => $flow_split_fastq,
      });
    push(@analyses, {
           -logic_name => 'bwa',
            -module        => 'ReseqTrack::HiveProcess::BWA',
            -parameters    => {
                program_file => $self->o('bwa_exe'),
                samtools => $self->o('samtools_exe'),
                reference => $self->o('reference'),
                branch_parameters_in => {
                    fastq => {key => 'FASTQ', inactivate => 1},
                    run_id => {key => 'RUN_ID', ascend => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -rc_name => '5Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ 'sort_chunks'],
            },
      });
    push(@analyses, {
            -logic_name => 'sort_chunks',
            -module        => 'ReseqTrack::HiveProcess::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                command => 'sort',
                create_index => 1,
                jvm_args => '-Xmx2g',
                branch_parameters_in => {
                    bam => {key => 'BAM', inactivate => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                  bai => 'BAI',
                },
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    if ($process_run_level) {
      push(@analyses, {
            -logic_name => 'decide_merge_chunks',
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
            -flow_into => $flow_decide_merge_chunks,
        });
      push(@analyses, {
            -logic_name => 'merge_chunks',
            -module        => 'ReseqTrack::HiveProcess::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                jvm_args => '-Xmx2g',
                command => 'merge',
                create_index => 1,
                branch_parameters_in => {
                    bam => {key => 'BAM', descend => 1, inactivate => 1},
                    bai => {key => 'BAI', descend => 1, inactivate => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                  bai => 'BAI',
                },
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    }
    push(@analyses, {
            -logic_name => 'decide_merge_runs',
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
                '2->A' => [ 'merge_runs'],
                'A->1' => [ 'mark_duplicates' ],
            },
      });
    push(@analyses, {
            -logic_name => 'merge_runs',
            -module        => 'ReseqTrack::HiveProcess::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                jvm_args => '-Xmx2g',
                command => 'merge',
                branch_parameters_in => {
                    bam => {key => 'BAM', descend => 1, inactivate => 1},
                    bai => {key => 'BAI', descend => 1, inactivate => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    push(@analyses, {
            -logic_name => 'decide_merge_libraries',
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
            -flow_into => $flow_decide_merge_libraries,
      });
    push(@analyses, {
            -logic_name => 'merge_libraries',
            -module        => 'ReseqTrack::HiveProcess::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                jvm_args => '-Xmx2g',
                command => 'merge',
                create_index => 1,
                branch_parameters_in => {
                    bam => {key => 'BAM', descend => 1, inactivate => 1},
                    bai => {key => 'BAI', descend => 1, inactivate => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                  bai => 'BAI',
                },
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    push(@analyses, {
            -logic_name => 'mark_duplicates',
            -module        => 'ReseqTrack::HiveProcess::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                jvm_args => '-Xmx4g',
                command => 'mark_duplicates',
                create_index => 1,
                branch_parameters_in => {
                    bam => {key => 'BAM', inactivate => 1, descend => 1},
                    bai => {key => 'BAI', inactivate => 1, descend => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                  bai => 'BAI',
                },
            },
            -rc_name => '5Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    push(@analyses, {
            -logic_name => 'realign',
            -module        => 'ReseqTrack::HiveProcess::RunBamImprovement',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'realign',
                reference => $self->o('reference'),
                gatk_dir => $self->o('gatk_dir'),
                known_sites_vcf => $self->o('known_indels_vcf'),
                gatk_module_options => {knowns_only => $realign_knowns_only},
                branch_parameters_in => {
                    bam => {key => 'BAM', inactivate => 1, descend => 1},
                    bai => {key => 'BAI', inactivate => 1, descend => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -rc_name => '5Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => $flow_realign,
      });
    if ($do_index_for_recalibrate) {
      push(@analyses, {
            -logic_name => 'index_for_recalibrate',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'index',
                branch_parameters_in => {
                    bam => {key => 'BAM', descend => 1},
                },
                branch_parameters_out => {
                  bai => 'BAI',
                },
            },
            -flow_into => {1 => ['recalibrate']},
            -rc_name => '200Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    }
    if ($recalibrate_level) {
      push(@analyses, {
            -logic_name => 'recalibrate',
            -module        => 'ReseqTrack::HiveProcess::RunBamImprovement',
            -meadow_type=> 'LOCAL',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'recalibrate',
                reference => $self->o('reference'),
                gatk_dir => $self->o('gatk_dir'),
                jvm_args => '-Xmx2g',
                known_sites_vcf => $self->o('known_snps_vcf'),
                branch_parameters_in => {
                    bam => {key => 'BAM', inactivate => 1, descend => 1},
                    bai => {key => 'BAI', inactivate => 1, descend => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => $flow_recalibrate,
      });
    }
    push(@analyses, {
            -logic_name => 'calmd',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'calmd',
                reference => $self->o('reference'),
                samtools_options => {input_sort_status => 'c'},
                branch_parameters_in => {
                    bam => {key => 'BAM', inactivate => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ 'tag_strip'],
            },
      });
    push(@analyses, {
            -logic_name => 'tag_strip',
            -module        => 'ReseqTrack::HiveProcess::RunSqueezeBam',
            -parameters => {
                'program_file' => $self->o('squeeze_exe'),
                'rm_OQ_fields' => 1,
                'rm_tag_types' => ['XM:i', 'XG:i', 'XO:i'],
                'branch_parameters_in' => {
                    bam => {key => 'BAM', inactivate => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -rc_name => '1Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => $flow_tag_strip,
      });

    push(@analyses, {
            -logic_name => 'reheader',
            -module        => 'ReseqTrack::HiveProcess::ReheaderBam',
            -parameters => {
                'samtools' => $self->o('samtools_exe'),
                'header_lines_file' => $self->o('header_lines_file'),
                'SQ_assembly' => $self->o('ref_assembly'),
                'SQ_species' => $self->o('ref_species'),
                'SQ_uri' => $self->o('reference_uri'),
                'branch_parameters_in' => {
                    bam => {key => 'BAM', inactivate => 1},
                    fastq => {key => 'SOURCE_FASTQ', descend => 1},
                },
                branch_parameters_out => {
                  bam => 'BAM',
                },
            },
            -rc_name => '1Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ 'rename'],
            },
      });
    push(@analyses, {
            -logic_name => 'rename',
            -module        => 'ReseqTrack::HiveProcess::RenameFile',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                analysis_label => $self->o('final_label'),
                suffix => 'bam',
                branch_parameters_in => {
                    old_filename => {key => 'BAM', inactivate => 1},
                },
                branch_parameters_out => {
                  new_filename => 'BAM',
                },
            },
            -flow_into => {
                1 => [ 'final_index'],
            },
      });
    push(@analyses, {
            -logic_name => 'final_index',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'index',
                branch_parameters_in => {
                    bam => 'BAM'
                },
                branch_parameters_out => {
                  bai => 'BAI',
                },
            },
            -flow_into => {1 => ['validate']},
            -rc_name => '1Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });
    push(@analyses, {
            -logic_name => 'validate',
            -module        => 'ReseqTrack::HiveProcess::RunValidateBam',
            -parameters => {
                'program_file' => $self->o('validate_bam_exe'),
                'branch_parameters_in' => {
                    bam => 'BAM',
                },
                branch_parameters_out => {
                  bas => 'BAS',
                },
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });

    return \@analyses;
}

1;

