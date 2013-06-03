

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
        'realign_intervals_file' => '',

        'reference_uri' => $self->o('reference'),

        'final_label' => $self->o('pipeline_name'),

        'realign_knowns_only' => 0,
        'recalibrate_level' => 2,

    };
}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands that will create and set up the Hive database.
                  In addition to the standard creation of the database and populating it with Hive tables and procedures it also creates two pipeline-specific tables used by Runnables to communicate.

=cut

sub pipeline_create_commands {
    my ($self) = @_;

    my $sql_1 = '
    CREATE TABLE reseqtrack_file (
      file_id int(10) unsigned NOT NULL AUTO_INCREMENT,
      name VARCHAR(64000) NOT NULL,
      PRIMARY KEY (file_id)
    )';

    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation

        $self->db_execute_command('pipeline_db', $sql_1),

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

        ##################  WRITE SOME CODE TO HANDLE THIS!!!! Do conversion on way in and way out.
        'file_param_names' => ['fastq', 'bam'],

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

    my @analyses;
    push(@analyses, {
            -logic_name    => 'studies_factory',
            -module        => 'ReseqTrack::HiveProcess::JobFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -input_ids => [{data_value => [split(',', $self->o('study_id'))]}],
            -flow_into => {
                2 => { 'samples_factory' => {'study_id' => '#data_value#'} },   # will create a semaphored fan of jobs
            },
      });
    push(@analyses, {
            -logic_name    => 'samples_factory',
            -module        => 'ReseqTrack::HiveProcess::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                type_branch => 'sample',
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
            },
            -flow_into => {
                '2->A' => [ 'find_source_fastqs' ],   # will create a semaphored fan of jobs
                'A->1' => [ 'decide_merge_runs'  ],   # will create a semaphored funnel job to wait for the fan to complete
            },
      });
    push(@analyses, {
            -logic_name    => 'find_source_fastqs',
            -module        => 'ReseqTrack::HiveProcess::ImportCollection',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters    => {
                collection_type => $self->o('type_fastq'),
                collection_name_param => 'run_id',
                output_name_param => 'fastq',
            },
            -flow_into => {
                1 => [ 'split_fastq', ':////accu?fastq=[]' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'split_fastq',
            -module        => 'ReseqTrack::HiveProcess::SplitFastq',
            -parameters    => {
                program_file => $self->o('split_exe'),
                max_bases => $self->o('chunk_max_bases'),
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
              '2->A' => ['bwa'],
              'A->1' => ['decide_merge_chunks'],
            }
      });
    push(@analyses, {
           -logic_name => 'bwa',
            -module        => 'ReseqTrack::HiveProcess::BWA',
            -parameters    => {
                program_file => $self->o('bwa_exe'),
                samtools => $self->o('samtools_exe'),
                reference => $self->o('reference'),
                delete_param => 'fastq',
            },
            -rc_name => '5Gb', # Note the 'hardened' version of BWA may need 8Gb RAM or more
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => ['sort_chunks'],
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
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]']
            },
      });
    push(@analyses, {
          -logic_name => 'decide_merge_chunks',
          -module        => 'ReseqTrack::HiveProcess::FlowDecider',
          -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -parameters => {
              add_params_to_output => ['bam', 'bai'],
              file_param => 'bam',
              flows_if_no_files => undef,
              flows_if_one_file => [1,3],
              flows_if_multiple_files => [1,2,3],
              realign_knowns_only => $self->o('realign_knowns_only'),
              recalibrate_level => $self->o('recalibrate_level'),
              require_true => {
                  1 => '#expr($realign_knowns_only || $recalibrate_level==1)expr#',
                  2 => '#expr($realign_knowns_only || $recalibrate_level==1)expr#',
                  3 => '#expr(!$realign_knowns_only && $recalibrate_level!=1)expr#',
              }
          },
            -flow_into => {
                '2->A' => [ 'merge_chunks' ],
                'A->1' => [ 'decide_realign_run_level' ],
                3 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
            },
      });
    push(@analyses, {
          -logic_name => 'merge_chunks',
          -module        => 'ReseqTrack::HiveProcess::RunPicard',
          -parameters => {
              picard_dir => $self->o('picard_dir'),
              jvm_args => '-Xmx2g',
              command => 'merge',
              create_index => 1,
              delete_param => ['bam', 'bai'],
          },
          -rc_name => '2Gb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
    });

    push(@analyses, {
            -logic_name => 'decide_realign_run_level',
            -module        => 'ReseqTrack::HiveProcess::FlowDecider',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                require_true => {2 => $self->o('realign_knowns_only')},
                file_param => 'bam',
                flows_if_no_files => undef,
                flows_if_one_file => [2,1],
                flows_if_multiple_files => [2,1],
            },
            -flow_into => {
                '2->A' => [ 'realign_knowns_only'],
                'A->1' => [ 'decide_recalibrate_run_level' ],
            },
      });
    push(@analyses, {
            -logic_name => 'realign_knowns_only',
            -module        => 'ReseqTrack::HiveProcess::RunBamImprovement',
            -parameters => {
                command => 'realign',
                reference => $self->o('reference'),
                gatk_dir => $self->o('gatk_dir'),
                known_sites_vcf => $self->o('known_indels_vcf'),
                intervals_file => $self->o('realign_intervals_file'),
                gatk_module_options => {knowns_only => 1},
                delete_param => ['bam', 'bai'],
            },
            -rc_name => '5Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ 'calmd_run_level'],
            },
      });
    push(@analyses, {
            -logic_name => 'calmd_run_level',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'calmd',
                reference => $self->o('reference'),
                samtools_options => {input_sort_status => 'c'},
                delete_param => ['bam'],
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ 'tag_strip_run_level'],
            },
      });
    push(@analyses, {
            -logic_name => 'tag_strip_run_level',
            -module        => 'ReseqTrack::HiveProcess::RunSqueezeBam',
            -parameters => {
                program_file => $self->o('squeeze_exe'),
                'rm_OQ_fields' => 1,
                'rm_tag_types' => ['XM:i', 'XG:i', 'XO:i'],
                delete_param => ['bam'],
            },
            -rc_name => '1Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ ':////accu?bam=[]'],
            },
      });
    push(@analyses, {
            -logic_name => 'decide_recalibrate_run_level',
            -module        => 'ReseqTrack::HiveProcess::FlowDecider',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                require_true => {2 => '#expr($recalibrate_level==1)expr#'},
                file_param => 'bam',
                flows_if_no_files => undef,
                flows_if_one_file => [1,2],
                flows_if_multiple_files => [2,1],
            },
            -flow_into => {
                '2->A' => [ 'decide_index_recalibrate_run_level'],
                'A->1' => [ 'dummy_process_1' ],
            },
      });
    push(@analyses, {
          -logic_name => 'decide_index_recalibrate_run_level',
          -module        => 'ReseqTrack::HiveProcess::FlowDecider',
          -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -parameters => {
              file_param => 'bai',
              flows_if_no_files => [1,2],
              flows_if_one_file => 1,
          },
            -flow_into => {
                '2->A' => [ 'index_recalibrate_run_level' ],
                'A->1' => [ 'recalibrate_run_level' ],
            },
      });
    push(@analyses, {
          -logic_name => 'index_recalibrate_run_level',
          -module        => 'ReseqTrack::HiveProcess::RunSamtools',
          -parameters => {
              program_file => $self->o('samtools_exe'),
              command => 'index',
          },
          -rc_name => '200Mb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
            -flow_into => {
                1 => [':////accu?bai=[]'],
            },
    });
    push(@analyses, {
          -logic_name => 'recalibrate_run_level',
          -module        => 'ReseqTrack::HiveProcess::RunBamImprovement',
          -parameters => {
              command => 'recalibrate',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              jvm_args => '-Xmx2g',
              known_sites_vcf => $self->o('known_snps_vcf'),
              delete_param => ['bam', 'bai'],
          },
          -rc_name => '2Gb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
    });
    push(@analyses, {
          -logic_name => 'dummy_process_1',
          -module        => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
    });
    push(@analyses, {
          -logic_name => 'decide_merge_runs',
          -module        => 'ReseqTrack::HiveProcess::FlowDecider',
          -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -parameters => {
              add_params_to_output => ['bam', 'bai', 'fastq'],
              file_param => 'bam',
              flows_if_no_files => undef,
              flows_if_one_file => 1,
              flows_if_multiple_files => [2,1],
          },
            -flow_into => {
                '2->A' => [ 'merge_runs' ],
                'A->1' => [ 'mark_duplicates' ],
                1 => [':////accu?fastq=[]'],
            },
    });
    push(@analyses, {
          -logic_name => 'merge_runs',
          -module        => 'ReseqTrack::HiveProcess::RunPicard',
          -parameters => {
              picard_dir => $self->o('picard_dir'),
              jvm_args => '-Xmx2g',
              command => 'merge',
              create_index => 1,
              delete_param => ['bam', 'bai'],
          },
          -rc_name => '2Gb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
    });
    push(@analyses, {
            -logic_name => 'mark_duplicates',
            -module        => 'ReseqTrack::HiveProcess::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                jvm_args => '-Xmx4g',
                command => 'mark_duplicates',
                create_index => 1,
                delete_param => ['bam', 'bai'],
            },
            -rc_name => '5Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
            },
    });
    push(@analyses, {
          -logic_name => 'decide_merge_libraries',
          -module        => 'ReseqTrack::HiveProcess::FlowDecider',
          -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -parameters => {
              add_params_to_output => ['bam', 'bai', 'fastq'],
              file_param => 'bam',
              flows_if_no_files => undef,
              flows_if_one_file => 1,
              flows_if_multiple_files => [2,1],
          },
            -flow_into => {
                '2->A' => [ 'merge_libraries' ],
                'A->1' => [ 'decide_realign_sample_level' ],
            },
    });
    push(@analyses, {
          -logic_name => 'merge_libraries',
          -module        => 'ReseqTrack::HiveProcess::RunPicard',
          -parameters => {
              picard_dir => $self->o('picard_dir'),
              jvm_args => '-Xmx2g',
              command => 'merge',
              create_index => 1,
              delete_param => ['bam', 'bai'],
          },
          -rc_name => '2Gb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
    });
    push(@analyses, {
            -logic_name => 'decide_realign_sample_level',
            -module        => 'ReseqTrack::HiveProcess::FlowDecider',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                realign_knowns_only => $self->o('realign_knowns_only'),
                require_true => {2 => '#expr(!$realign_knowns_only)expr#'},
                file_param => 'bam',
                flows_if_no_files => undef,
                flows_if_one_file => [2,1],
                flows_if_multiple_files => [2,1],
            },
            -flow_into => {
                '2->A' => [ 'realign_full'],
                'A->1' => [ 'decide_recalibrate_sample_level' ],
            },
      });
    push(@analyses, {
            -logic_name => 'realign_full',
            -module        => 'ReseqTrack::HiveProcess::RunBamImprovement',
            -parameters => {
                command => 'realign',
                reference => $self->o('reference'),
                gatk_dir => $self->o('gatk_dir'),
                known_sites_vcf => $self->o('known_indels_vcf'),
                gatk_module_options => {knowns_only => 0},
                delete_param => ['bam', 'bai'],
            },
            -rc_name => '5Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ 'calmd_sample_level'],
            },
      });
    push(@analyses, {
            -logic_name => 'calmd_sample_level',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'calmd',
                reference => $self->o('reference'),
                samtools_options => {input_sort_status => 'c'},
                delete_param => ['bam'],
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ 'tag_strip_sample_level'],
            },
      });
    push(@analyses, {
            -logic_name => 'tag_strip_sample_level',
            -module        => 'ReseqTrack::HiveProcess::RunSqueezeBam',
            -parameters => {
                program_file => $self->o('squeeze_exe'),
                'rm_OQ_fields' => 1,
                'rm_tag_types' => ['XM:i', 'XG:i', 'XO:i'],
                delete_param => ['bam'],
            },
            -rc_name => '1Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
            },
      });
    push(@analyses, {
            -logic_name => 'decide_recalibrate_sample_level',
            -module        => 'ReseqTrack::HiveProcess::FlowDecider',
            -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                recalibrate_level => $self->o('recalibrate_level'),
                require_true => {2 => '#expr($recalibrate_level==2)expr#'},
                file_param => 'bam',
                flows_if_no_files => undef,
                flows_if_one_file => [1,2],
                flows_if_multiple_files => [1,2],
            },
            -flow_into => {
                '2->A' => [ 'decide_index_recalibrate_sample_level'],
                'A->1' => [ 'reheader' ],
            },
      });
    push(@analyses, {
          -logic_name => 'decide_index_recalibrate_sample_level',
          -module        => 'ReseqTrack::HiveProcess::FlowDecider',
          -meadow_type=> 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
          -parameters => {
              file_param => 'bai',
              flows_if_no_files => [1,2],
              flows_if_one_file => 1,
          },
            -flow_into => {
                '2->A' => [ 'index_recalibrate_sample_level' ],
                'A->1' => [ 'recalibrate_sample_level' ],
            },
      });
    push(@analyses, {
          -logic_name => 'index_recalibrate_sample_level',
          -module        => 'ReseqTrack::HiveProcess::RunSamtools',
          -parameters => {
              program_file => $self->o('samtools_exe'),
              command => 'index',
          },
          -rc_name => '200Mb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
            -flow_into => {
                1 => [':////accu?bai=[]'],
            },
    });
    push(@analyses, {
          -logic_name => 'recalibrate_sample_level',
          -module        => 'ReseqTrack::HiveProcess::RunBamImprovement',
          -parameters => {
              command => 'recalibrate',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              jvm_args => '-Xmx2g',
              known_sites_vcf => $self->o('known_snps_vcf'),
              delete_param => ['bam', 'bai'],
          },
          -rc_name => '2Gb',
          -analysis_capacity  =>  50,  # use per-analysis limiter
          -hive_capacity  =>  -1,
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
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
                delete_param => ['bam', 'bai'],
            },
            -rc_name => '1Gb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
            -flow_into => {
                1 => {'rename' => {'file' => '#bam#'} },
            },
      });
    push(@analyses, {
            -logic_name => 'rename',
            -module        => 'ReseqTrack::HiveProcess::RenameFile',
            -meadow_type => 'LOCAL',     # do not bother the farm with such a simple task (and get it done faster)
            -parameters => {
                analysis_label => $self->o('final_label'),
                suffix => 'bam',
            },
            -flow_into => {
                1 => {'final_index' => {'bam' => '#file#'} },
            },
      });
    push(@analyses, {
            -logic_name => 'final_index',
            -module        => 'ReseqTrack::HiveProcess::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'index',
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
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  50,  # use per-analysis limiter
            -hive_capacity  =>  -1,
      });


    return \@analyses;
}

1;

