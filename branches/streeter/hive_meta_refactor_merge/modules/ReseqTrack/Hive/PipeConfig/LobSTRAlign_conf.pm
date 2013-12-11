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
      -sample_label, (default sample_source_id), which column from the sample table should be used to name output files and directories
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


package ReseqTrack::Hive::PipeConfig::LobSTRAlign_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');


sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },

        'pipeline_name' => 'align',

        seeding_module => 'ReseqTrack::Hive::PipeSeed::BasePipeSeed',
        seeding_options => {
            output_columns => $self->o('sample_columns'),
            output_attributes => $self->o('sample_attributes'),
            require_columns => $self->o('require_sample_columns'),
            exclude_columns => $self->o('exclude_sample_columns'),
            require_attributes => $self->o('require_sample_attributes'),
            exclude_attributes => $self->o('exclude_sample_attributes'),
          },

        'chunk_max_reads'    => 5000000,
        'type_fastq'    => 'FILTERED_FASTQ',
        'split_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/split/split',
        'lobstr_exe' => '/nfs/1000g-work/G1K/work/bin/lobstr/bin/lobSTR',
        'samtools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/samtools',
        'picard_dir' => '/nfs/1000g-work/G1K/work/bin/picard',

        #various options for overriding defaults in reheadering bam
        'dict_file' => undef,
        'reference_uri' => undef,
        'ref_assembly' => undef,
        'ref_species' => undef,
        'header_lines_file' => undef,

        'bam_type' => undef,
        'bai_type' => undef,

        'lobstr_options' => { threads => 4 },

        'RGSM' => '#sample_source_id#',
        'RGPU' => '#run_source_id#',

        'sample_attributes' => [],
        'sample_columns' => ['sample_id', 'sample_source_id', 'sample_alias'],
        'run_attributes' => [],
        'run_columns' => ['run_source_id', 'center_name', 'run_alias'],
        'study_attributes' => [],
        'study_columns' => ['study_source_id'],
        'experiment_attributes' => [],
        'experiment_columns' => ['instrument_platform', 'paired_nominal_length'],

        require_run_attributes => {},
        require_experiment_attributes => {},
        require_study_attributes => {},
        require_sample_attributes => {},
        exclude_run_attributes => {},
        exclude_experiment_attributes => {},
        exclude_study_attributes => {},
        exclude_sample_attributes => {},
        require_experiment_columns => { instrument_platform => ['ILLUMINA'], },
        require_run_columns => { status => ['public', 'private'], },
        require_study_columns => {},
        require_sample_columns => {},
        exclude_sample_columns => {},

        final_output_dir => $self->o('root_output_dir'),
        final_output_layout => '#sample_source_id#/alignment',
        name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method => 'basic',
        name_file_params => {
            new_dir => '#final_output_dir#/#final_output_layout#',
            new_basename => '#sample_source_id#.lobstr',
            add_datestamp => 1,
            suffix => ['.bam', '.bam.bai'],
          },

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
        dir_label_params => ["sample_source_id", "library_name", "run_source_id", "chunk_label"],
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
            '6Gb' => { 'LSF' => '-C0 -M6000 -q production -R"select[mem>6000] rusage[mem=6000]"' },
            '6Gb_4t' => { 'LSF' => '-C0 -M6000 -q production -R"select[mem>6000] rusage[mem=6000]" -n4' },
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
                seeding_module => $self->o('seeding_module'),
                seeding_options => $self->o('seeding_options'),
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
                require_experiment_columns => $self->o('require_experiment_columns'),
                require_study_columns => $self->o('require_study_columns'),
                require_experiment_attributes => $self->o('require_experiment_attributes'),
                require_study_attributes => $self->o('require_study_attributes'),
                exclude_experiment_attributes => $self->o('exclude_experiment_attributes'),
                exclude_study_attributes => $self->o('exclude_study_attributes'),
            },
            -flow_into => {
                '2->A' => [ 'runs_factory' ],
                'A->1' => [ 'decide_merge_bams'  ],
            },
      });
    push(@analyses, {
            -logic_name    => 'runs_factory',
            -module        => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                factory_type => 'run',
                require_experiment_columns => $self->o('require_experiment_columns'),
                require_study_columns => $self->o('require_study_columns'),
                require_run_columns => $self->o('require_run_columns'),
                require_experiment_attributes => $self->o('require_experiment_attributes'),
                require_study_attributes => $self->o('require_study_attributes'),
                require_run_attributes => $self->o('require_run_attributes'),
                exclude_experiment_attributes => $self->o('exclude_experiment_attributes'),
                exclude_study_attributes => $self->o('exclude_study_attributes'),
                exclude_run_attributes => $self->o('exclude_run_attributes'),

                output_run_columns => $self->o('run_columns'),
                output_study_columns => $self->o('study_columns'),
                output_experiment_columns => $self->o('experiment_columns'),
                output_run_attributes => $self->o('run_attributes'),
                output_study_attributes => $self->o('study_attributes'),
                output_experiment_attributes => $self->o('experiment_attributes'),
            },
            -flow_into => {
                '2->A' => [ 'find_source_fastqs' ],
                'A->1' => [ 'collect_bams2'  ],
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
                reseqtrack_options => {
                  flows_do_count_param => 'fastq',
                  flows_do_count => { 1 => '1+', },
                },
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
              '2->A' => {'lobstr' => {'chunk_label' => '#expr($run_source_id.".".$chunk)expr#', 'fastq' => '#fastq#'}},
              'A->1' => ['collect_bams1'],
            }
      });
    push(@analyses, {
           -logic_name => 'lobstr',
            -module        => 'ReseqTrack::Hive::Process::LobSTR',
            -parameters    => {
                program_file => $self->o('lobstr_exe'),
                ref_index_prefix => $self->o('lobstr_ref_index_prefix'),
                options => $self->o('lobstr_options'),
                reseqtrack_options => {
                  delete_param => 'fastq',
                },
              },
            -rc_name => '6Gb_4t',
            -hive_capacity  =>  100,
            -flow_into => {
                1 => ['sort_chunks'],
            },
      });
    push(@analyses, {
            -logic_name => 'sort_chunks',
            -module        => 'ReseqTrack::Hive::Process::RunSamtools',
            -parameters => {
                program_file => $self->o('samtools_exe'),
                command => 'sort',
                reseqtrack_options => {
                  delete_param => 'bam',
                },
            },
            -rc_name => '2Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [ 'fix_read_group']
            },
      });
    push(@analyses, {
            -logic_name => 'fix_read_group',
            -module        => 'ReseqTrack::Hive::Process::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                command => 'add_or_replace_read_groups',
                create_index => 0,
                jvm_args => '-Xmx2g',
                options => {read_group_fields => {
                  ID => '#run_source_id#',
                  LB => '#library_name#',
                  CN => '#center_name#',
                  PI => '#paired_nominal_length#',
                  SM => $self->o('RGSM'),
                  DS => '#study_source_id#',
                  PU => $self->o('RGPU'),
                  PL => '#instrument_platform#',
                }},
                reseqtrack_options => {
                  delete_param => 'bam',
                },
            },
            -rc_name => '2Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [ ':////accu?bam=[]']
            },
      });
    push(@analyses, {
          -logic_name => 'collect_bams1',
          -module        => 'ReseqTrack::Hive::Process::BaseProcess',
          -meadow_type=> 'LOCAL',
          -parameters => {
          },
          -flow_into => {
                1 => [ ':////accu?bam=[]'],
          },
      });
    push(@analyses, {
          -logic_name => 'collect_bams2',
          -module        => 'ReseqTrack::Hive::Process::BaseProcess',
          -meadow_type=> 'LOCAL',
          -parameters => {
          },
          -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?fastq=[]'],
          },
      });
    push(@analyses, {
          -logic_name => 'decide_merge_bams',
          -module        => 'ReseqTrack::Hive::Process::BaseProcess',
          -meadow_type=> 'LOCAL',
          -parameters => {
              reseqtrack_options => {
                denestify => 'bam',
                flows_non_factory => [1,2,7],
                flows_do_count_param => 'bam',
                flows_do_count => {
                          1 => '1+',
                          2 => '2+',
                          7 => '0',
                        },
              },
          },
            -flow_into => {
                '2->A' => [ 'merge_bams' ],
                'A->1' => [ 'reheader' ],
                '7' => [ 'mark_seed_futile' ],
            },
    });
    push(@analyses, {
          -logic_name => 'merge_bams',
          -module        => 'ReseqTrack::Hive::Process::RunPicard',
          -parameters => {
              picard_dir => $self->o('picard_dir'),
              jvm_args => '-Xmx2g',
              command => 'merge',
              create_index => 0,
              reseqtrack_options => {
                denestify => 'bam',
                delete_param => 'bam',
              },
          },
          -rc_name => '2Gb',
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
                reseqtrack_options => {
                  denestify => 'bam',
                  delete_param => 'bam',
                },
            },
            -rc_name => '200Mb',
            -hive_capacity  =>  200,
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
            -flow_into => {1 => ['store_bam']},
            -rc_name => '200Mb',
            -hive_capacity  =>  200,
      });
    push(@analyses, {
            -logic_name    => 'store_bam',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -parameters    => {
              type => $self->o('bam_type'),
              file => '#bam#',
              name_file_module => $self->o('name_file_module'),
              name_file_method => $self->o('name_file_method'),
              name_file_params => $self->o('name_file_params'),
              final_output_dir => $self->o('final_output_dir'),
              final_output_layout => $self->o('final_output_layout'),
            },
            -rc_name => '200Mb',
            -hive_capacity  =>  200,
            -flow_into => {1 => {'store_bai' => {'bam' => '#file#'}}},
      });
    push(@analyses, {
            -logic_name    => 'store_bai',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -parameters    => {
              type => $self->o('bai_type'),
              file => '#bai#',
              name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
              name_file_method => 'basic',
              name_file_params => {new_full_path => '#bam#.bai'},
            },
            -rc_name => '200Mb',
            -hive_capacity  =>  200,
            -flow_into => {1 => {'mark_seed_complete' => {'bai' => '#file#'}}},
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

