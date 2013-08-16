=head1 NAME

 ReseqTrack::Hive::PipeConfig::Alignment_conf

=head1 SYNOPSIS

  This is a pipeline for aligning fastq reads to make a bam file.

  Pipeline MUST be seeded by the sample table of a ReseqTrack database.
  Bam/bai/bas files will be created for each sample and stored in the ReseqTrack database

  Here is an example pipeline configuration to load using reseqtrack/scripts/pipeline/load_pipeline_from_conf.pl

[alignment]
table_name=sample
seeding_module=ReseqTrack::Hive::PipeSeed::Default
config_module=ReseqTrack::Hive::PipeConfig::Alignment_conf
config_options=-root_output_dir /path/to/dir
config_options=-reference /path/to/human.fa
config_options=-known_indels_vcf /nfs/1000g-archive/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz
config_options=-known_snps_vcf /nfs/1000g-archive/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz
config_options=-realign_intervals_file /path/to/realigner.intervals
config_options=-realign_knowns_only 1
config_options=-recalibrate_level 1
config_options=-allowed_strategy WGS

  Options that MUST be specified in the pipeline.config_options table/column of your ReseqTrack database:

      -reference, fasta file of your reference genome.  Should be indexed for bwa and should have a .fai and .dict
      -allowed_strategy, e.g. WGS or WXS, taken from the experiment table.

  Options that have defaults but you will often want to set them in your pipeline.cofig_options table/column:

      -root_output_dir, (default is your current directory)
      -type_fastq, type of fastq files to look for in the reseqtrack database, default FILTERED_FASTQ
      -final_label, used to name your final output files (default is your pipeline name)
      -sample_label, (default source_id), which column from the sample table should be used to name output files and directories
      -sample_group_attribute, (default POPULATION), output bams will be grouped in directories using this attribute

      -chunk_max_reads, (default 5000000) controls how fastq files are split up into chunks for parallel alignment

      Recalibration and Realignment options:
      -realign_knowns_only, boolean, default 0.  You can choose between full indel realignment (1) or faster realignment around known indels only (0)
      -recalibrate_level, can be 0 (don't recalibrate), 1 (fast recalibration at lane level, e.g. 1000genomes), 2 (slower recalibration at sample level for better accuracy)
      -known_indels_vcf, used for indel realignment (default undefined).  Optional if realign_knowns_only=0; mandatory if realign_knowns_only=1.
      -known_snps_vcf, used for recalibration (default undefined). Mandatory if recalibrate_level != 0.
      -realign_intervals_file, should be given if realign_known_only=1.  Can be generated using gatk RealignerTargetCreator

      Various options for reheadering a bam file:
      -header_lines_file should contain any @PG and @CO lines you want written to your bam header.
      -dict_file.  @SQ lines from this file will be written to the bam header.  Default is to use the dict file associated with your reference file
      -reference_uri. (undef) Used to override default in the @SQ lines.
      -ref_species. (undef) Used to override default in the @SQ lines.
      -ref_assembly. (undef) Used to override default in the @SQ lines.

      Paths of executables:
      -split_exe, (default is to work it out from your environment variable $RESEQTRACK)
      -validate_bam_exe, (default is to work it out from your environment variable $RESEQTRACK)
      -bwa_exe, (default /nfs/1000g-work/G1K/work/bin/bwa/bwa)
      -samtools_exe => (default /nfs/1000g-work/G1K/work/bin/samtools/samtools)
      -squeeze_exe => (default /nfs/1000g-work/G1K/work/bin/bamUtil/bin/bam)
      -gatk_dir => (default /nfs/1000g-work/G1K/work/bin/gatk/dist/)
      -picard_dir => (default /nfs/1000g-work/G1K/work/bin/picard)

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


package ReseqTrack::Hive::PipeConfig::Alignment_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');


sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },

        'pipeline_name' => 'align',

        'chunk_max_reads'    => 5000000,
        'type_fastq'    => 'FILTERED_FASTQ',
        'split_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/split/split',
        'validate_bam_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/validate_bam/validate_bam',
        'bwa_exe' => '/nfs/1000g-work/G1K/work/bin/bwa/bwa',
        'samtools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/samtools',
        'squeeze_exe' => '/nfs/1000g-work/G1K/work/bin/bamUtil/bin/bam',
        'gatk_dir' => '/nfs/1000g-work/G1K/work/bin/gatk/dist/',
        'picard_dir' => '/nfs/1000g-work/G1K/work/bin/picard',
        'known_indels_vcf' => undef,
        'known_snps_vcf' => undef,
        'realign_intervals_file' => undef,


        #various options for overriding defaults in reheadering bam
        'dict_file' => undef,
        'reference_uri' => undef,
        'ref_assembly' => undef,
        'ref_species' => undef,
        'header_lines_file' => undef,

        'final_label' => $self->o('pipeline_name'),

        'sample_label' => 'source_id',
        'sample_group_attribute' => 'POPULATION',

        'realign_knowns_only' => 0,
        'recalibrate_level' => 2,

        'bam_type' => undef,
        'bai_type' => undef,
        'bas_type' => undef,

    };
}


sub pipeline_create_commands {
    my ($self) = @_;

    return [
        @{$self->SUPER::pipeline_create_commands},
    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},
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
            '6Gb' => { 'LSF' => '-C0 -M6000 -q production -R"select[mem>6000] rusage[mem=6000]"' },
            '8Gb' => { 'LSF' => '-C0 -M8000 -q production -R"select[mem>8000] rusage[mem=8000]"' },
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
                sample_label => $self->o('sample_label'),
                sample_group_attribute => $self->o('sample_group_attribute'),
                seed_label => ['#sample_group_attribute#', '#sample_label#'],
                output_columns => ['sample_id', '#sample_label#'],
                output_attributes => '#sample_group_attribute#',
                temp_param_sub => { 2 => [['#sample_label#','undef'], ['#sample_group_attribute#', 'undef'], ['seed_time', 'undef']]}, # temporary hack pending updates to hive code
            },
            -flow_into => {
                2 => [ 'libraries_factory' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'libraries_factory',
            -module        => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                factory_type => 'library',
                allowed_strategy => $self->o('allowed_strategy'),
            },
            -flow_into => {
                '2->A' => [ 'runs_factory' ],
                'A->1' => [ 'decide_merge_libraries'  ],
            },
      });
    push(@analyses, {
            -logic_name    => 'runs_factory',
            -module        => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                factory_type => 'run',
                temp_param_sub => { 1 => [['run_source_id','undef']]}, # temporary hack pending updates to hive code
            },
            -flow_into => {
                '2->A' => [ 'find_source_fastqs' ],
                'A->1' => [ 'decide_mark_duplicates'  ],
            },
      });
    push(@analyses, {
            -logic_name    => 'find_source_fastqs',
            -module        => 'ReseqTrack::Hive::Process::ImportCollection',
            -meadow_type => 'LOCAL',
            -parameters    => {
                collection_type => $self->o('type_fastq'),
                collection_name => '#run_source_id#',
                output_param => 'fastq',
              temp_param_sub => { 1 => [['run_source_id','undef']]}, # temporary hack pending updates to hive code
            },
            -flow_into => {
                1 => [ 'split_fastq', ':////accu?fastq=[]' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'split_fastq',
            -module        => 'ReseqTrack::Hive::Process::SplitFastq',
            -parameters    => {
                program_file => $self->o('split_exe'),
                max_reads => $self->o('chunk_max_reads'),
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -flow_into => {
              '2->A' => ['bwa'],
              'A->1' => ['decide_merge_chunks'],
            }
      });
    push(@analyses, {
           -logic_name => 'bwa',
            -module        => 'ReseqTrack::Hive::Process::BWA',
            -parameters    => {
                program_file => $self->o('bwa_exe'),
                samtools => $self->o('samtools_exe'),
                reference => $self->o('reference'),
                delete_param => 'fastq',
            },
            -rc_name => '8Gb', # Note the 'hardened' version of BWA may need 8Gb RAM or more
            -hive_capacity  =>  100,
            -flow_into => {
                1 => ['sort_chunks'],
            },
      });
    push(@analyses, {
            -logic_name => 'sort_chunks',
            -module        => 'ReseqTrack::Hive::Process::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                command => 'sort',
                create_index => 1,
                jvm_args => '-Xmx2g',
                delete_param => 'bam',
            },
            -rc_name => '2Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]']
            },
      });
    push(@analyses, {
          -logic_name => 'decide_merge_chunks',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              realign_knowns_only => $self->o('realign_knowns_only'),
              recalibrate_level => $self->o('recalibrate_level'),
              files => '#bam#',
              require_file_count => {
                        1 => '1+',
                        2 => '2+',
                        3 => '1+',
                      },
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
          -module        => 'ReseqTrack::Hive::Process::RunPicard',
          -parameters => {
              picard_dir => $self->o('picard_dir'),
              jvm_args => '-Xmx2g',
              command => 'merge',
              create_index => 1,
              delete_param => ['bam', 'bai'],
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
    });

    push(@analyses, {
            -logic_name => 'decide_realign_run_level',
            -module        => 'ReseqTrack::Hive::Process::FlowDecider',
            -meadow_type=> 'LOCAL',
            -parameters => {
                require_true => {2 => $self->o('realign_knowns_only'), 1 => 1},
            },
            -flow_into => {
                '2->A' => [ 'realign_knowns_only'],
                'A->1' => [ 'decide_recalibrate_run_level' ],
            },
      });
    push(@analyses, {
            -logic_name => 'realign_knowns_only',
            -module        => 'ReseqTrack::Hive::Process::RunBamImprovement',
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
            -hive_capacity  =>  100,
            -flow_into => {
                1 => [ 'calmd_run_level'],
            },
      });
    push(@analyses, {
            -logic_name => 'calmd_run_level',
            -module        => 'ReseqTrack::Hive::Process::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'calmd',
                reference => $self->o('reference'),
                samtools_options => {input_sort_status => 'c'},
                delete_param => ['bam'],
            },
            -rc_name => '2Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
            },
      });
    push(@analyses, {
            -logic_name => 'decide_recalibrate_run_level',
            -module        => 'ReseqTrack::Hive::Process::FlowDecider',
            -meadow_type=> 'LOCAL',
            -parameters => {
                recalibrate_level => $self->o('recalibrate_level'),
                require_true => {2 => '#expr($recalibrate_level==1)expr#', 1 => 1},
            },
            -flow_into => {
                '2->A' => [ 'decide_index_recalibrate_run_level'],
                'A->1' => [ 'decide_tag_strip_run_level'],
            },
      });
    push(@analyses, {
          -logic_name => 'decide_index_recalibrate_run_level',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              count_files => 1,
              files => '#bai#',
              flows_if_no_files => 1,
              flows_if_one_file => [1,2],
              require_file_count => {
                        1 => '0+',
                        2 => '0',
                      },
          },
            -flow_into => {
                '2->A' => [ 'index_recalibrate_run_level' ],
                'A->1' => [ 'recalibrate_run_level' ],
            },
      });
    push(@analyses, {
          -logic_name => 'index_recalibrate_run_level',
          -module        => 'ReseqTrack::Hive::Process::RunSamtools',
          -parameters => {
              program_file => $self->o('samtools_exe'),
              command => 'index',
          },
          -rc_name => '200Mb',
          -hive_capacity  =>  200,
            -flow_into => {
                1 => [':////accu?bai=[]'],
            },
    });
    push(@analyses, {
          -logic_name => 'recalibrate_run_level',
          -module        => 'ReseqTrack::Hive::Process::RunBamImprovement',
          -parameters => {
              command => 'recalibrate',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              jvm_args => '-Xmx2g',
              known_sites_vcf => $self->o('known_snps_vcf'),
              delete_param => ['bam', 'bai'],
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
    });
    push(@analyses, {
            -logic_name => 'decide_tag_strip_run_level',
            -module        => 'ReseqTrack::Hive::Process::FlowDecider',
            -meadow_type=> 'LOCAL',
            -parameters => {
                recalibrate_level => $self->o('recalibrate_level'),
                realign_knowns_only => $self->o('realign_knowns_only'),
                require_true => {2 => '#expr($realign_knowns_only && $recalibrate_level<2)expr#',
                                 1 => '#expr(!$realign_knowns_only || $recalibrate_level==2)expr#'},
            },
            -flow_into => {
                2 => [ 'tag_strip_run_level'],
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
            },
      });
    push(@analyses, {
            -logic_name => 'tag_strip_run_level',
            -module        => 'ReseqTrack::Hive::Process::RunSqueezeBam',
            -parameters => {
                program_file => $self->o('squeeze_exe'),
                'rm_OQ_fields' => 1,
                'rm_tag_types' => ['XM:i', 'XG:i', 'XO:i'],
                delete_param => ['bam'],
            },
            -rc_name => '1Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
            },
      });
    push(@analyses, {
          -logic_name => 'decide_mark_duplicates',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              files => '#bam#',
              require_file_count => { 1 => '1+'},
          },
            -flow_into => {
                1 => [ 'mark_duplicates', ':////accu?fastq=[]'],
            },
      });
    push(@analyses, {
            -logic_name => 'mark_duplicates',
            -module        => 'ReseqTrack::Hive::Process::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                jvm_args => '-Xmx4g',
                command => 'mark_duplicates',
                create_index => 1,
                delete_param => ['bam', 'bai'],
            },
            -rc_name => '5Gb',
            -hive_capacity  =>  100,
            -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
            },
    });
    push(@analyses, {
          -logic_name => 'decide_merge_libraries',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              files => '#bam#',
              require_file_count => {
                        1 => '1+',
                        2 => '2+',
                        3 => '0',
                      }
          },
            -flow_into => {
                '2->A' => [ 'merge_libraries' ],
                'A->1' => [ 'decide_realign_sample_level' ],
                '3' => [ 'mark_seed_futile' ],
            },
    });
    push(@analyses, {
          -logic_name => 'merge_libraries',
          -module        => 'ReseqTrack::Hive::Process::RunPicard',
          -parameters => {
              picard_dir => $self->o('picard_dir'),
              jvm_args => '-Xmx2g',
              command => 'merge',
              create_index => 1,
              delete_param => ['bam', 'bai'],
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
    });
    push(@analyses, {
            -logic_name => 'decide_realign_sample_level',
            -module        => 'ReseqTrack::Hive::Process::FlowDecider',
            -meadow_type=> 'LOCAL',
            -parameters => {
                realign_knowns_only => $self->o('realign_knowns_only'),
                require_true => {2 => '#expr(!$realign_knowns_only)expr#', 1 => 1},
            },
            -flow_into => {
                '2->A' => [ 'realign_full'],
                'A->1' => [ 'decide_recalibrate_sample_level' ],
            },
      });
    push(@analyses, {
            -logic_name => 'realign_full',
            -module        => 'ReseqTrack::Hive::Process::RunBamImprovement',
            -parameters => {
                command => 'realign',
                reference => $self->o('reference'),
                gatk_dir => $self->o('gatk_dir'),
                known_sites_vcf => $self->o('known_indels_vcf'),
                gatk_module_options => {knowns_only => 0},
                delete_param => ['bam', 'bai'],
            },
            -rc_name => '5Gb',
            -hive_capacity  =>  100,
            -flow_into => {
                1 => [ 'calmd_sample_level'],
            },
      });
    push(@analyses, {
            -logic_name => 'calmd_sample_level',
            -module        => 'ReseqTrack::Hive::Process::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'calmd',
                reference => $self->o('reference'),
                samtools_options => {input_sort_status => 'c'},
                delete_param => ['bam'],
            },
            -rc_name => '2Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
            },
      });
    push(@analyses, {
            -logic_name => 'decide_recalibrate_sample_level',
            -module        => 'ReseqTrack::Hive::Process::FlowDecider',
            -meadow_type=> 'LOCAL',
            -parameters => {
                recalibrate_level => $self->o('recalibrate_level'),
                require_true => {2 => '#expr($recalibrate_level==2)expr#', 1=>1},
            },
            -flow_into => {
                '2->A' => [ 'decide_index_recalibrate_sample_level'],
                'A->1' => [ 'decide_tag_strip_sample_level' ],
            },
      });
    push(@analyses, {
          -logic_name => 'decide_index_recalibrate_sample_level',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              files => '#bai#',
              require_file_count => {
                        1 => '0+',
                        2 => '0',
                      },
          },
            -flow_into => {
                '2->A' => [ 'index_recalibrate_sample_level' ],
                'A->1' => [ 'recalibrate_sample_level' ],
            },
      });
    push(@analyses, {
          -logic_name => 'index_recalibrate_sample_level',
          -module        => 'ReseqTrack::Hive::Process::RunSamtools',
          -parameters => {
              program_file => $self->o('samtools_exe'),
              command => 'index',
          },
          -rc_name => '200Mb',
          -hive_capacity  =>  200,
            -flow_into => {
                1 => [':////accu?bai=[]'],
            },
    });
    push(@analyses, {
          -logic_name => 'recalibrate_sample_level',
          -module        => 'ReseqTrack::Hive::Process::RunBamImprovement',
          -parameters => {
              command => 'recalibrate',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              jvm_args => '-Xmx2g',
              known_sites_vcf => $self->o('known_snps_vcf'),
              delete_param => ['bam', 'bai'],
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?bam=[]', ':////accu?bai=[]'],
          },
    });
    push(@analyses, {
            -logic_name => 'decide_tag_strip_sample_level',
            -module        => 'ReseqTrack::Hive::Process::FlowDecider',
            -meadow_type=> 'LOCAL',
            -parameters => {
                recalibrate_level => $self->o('recalibrate_level'),
                realign_knowns_only => $self->o('realign_knowns_only'),
                require_true => {2 => '#expr(!$realign_knowns_only || $recalibrate_level==2)expr#', 1 => 1},
            },
            -flow_into => {
                '2->A' => [ 'tag_strip_sample_level'],
                'A->1' => [ 'reheader'],
            },
      });
    push(@analyses, {
            -logic_name => 'tag_strip_sample_level',
            -module        => 'ReseqTrack::Hive::Process::RunSqueezeBam',
            -parameters => {
                program_file => $self->o('squeeze_exe'),
                'rm_OQ_fields' => 1,
                'rm_tag_types' => ['XM:i', 'XG:i', 'XO:i'],
                delete_param => ['bam'],
            },
            -rc_name => '1Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [ ':////accu?bam=[]'],
            },
      });
    push(@analyses, {
            -logic_name => 'reheader',
            -module        => 'ReseqTrack::Hive::Process::ReheaderBam',
            -parameters => {
                'samtools' => $self->o('samtools_exe'),
                'header_lines_file' => $self->o('header_lines_file'),
                'dict_file' => $self->o('dict_file'),
                'reference' => $self->o('reference'),
                'SQ_assembly' => $self->o('ref_assembly'),
                'SQ_species' => $self->o('ref_species'),
                'SQ_uri' => $self->o('reference_uri'),
                delete_param => ['bam', 'bai'],
            },
            -rc_name => '1Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => ['rename'],
            },
      });
    push(@analyses, {
            -logic_name => 'rename',
            -module        => 'ReseqTrack::Hive::Process::RenameFile',
            -meadow_type => 'LOCAL',
            -parameters => {
                analysis_label => $self->o('final_label'),
                suffix => 'bam',
                file_param_name => 'bam',
            },
            -flow_into => {
                1 => ['final_index'],
            },
      });
    push(@analyses, {
            -logic_name => 'final_index',
            -module        => 'ReseqTrack::Hive::Process::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'index',
            },
            -flow_into => {1 => ['validate']},
            -rc_name => '1Gb',
            -hive_capacity  =>  200,
      });
    push(@analyses, {
            -logic_name => 'validate',
            -module        => 'ReseqTrack::Hive::Process::RunValidateBam',
            -parameters => {
                'program_file' => $self->o('validate_bam_exe'),
            },
            -rc_name => '200Mb',
            -hive_capacity  =>  200,
            -flow_into => {1 => ['store_bam']},
      });
    push(@analyses, {
            -logic_name    => 'store_bam',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('bam_type'),
              file => '#bam#',
            },
            -flow_into => {1 => ['store_bas']},
      });
    push(@analyses, {
            -logic_name    => 'store_bas',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('bas_type'),
              file => '#bas#',
            },
            -flow_into => {1 => ['store_bai']},
      });
    push(@analyses, {
            -logic_name    => 'store_bai',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('bai_type'),
              file => '#bai#',
            },
            -flow_into => {1 => ['mark_seed_complete']},
      });
    push(@analyses, {
            -logic_name    => 'mark_seed_complete',
            -module        => 'ReseqTrack::Hive::Process::UpdateSeed',
            -parameters    => {
              is_complete  => 1,
            },
            -meadow_type => 'LOCAL',
      });
    push(@analyses, {
            -logic_name    => 'mark_seed_futile',
            -module        => 'ReseqTrack::Hive::Process::UpdateSeed',
            -parameters    => {
              is_futile  => 1,
            },
            -meadow_type => 'LOCAL',
      });


    return \@analyses;
}

1;

