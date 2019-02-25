=head1 NAME

 ReseqTrack::Hive::PipeConfig::Alignment_conf

=head1 SYNOPSIS

  This is a pipeline for aligning fastq reads to make a bam and a cram file.

  Pipeline must be seeded by the sample table of a ReseqTrack database.
  Bam/bai/bas files will be created for each sample and stored in the ReseqTrack database

  Here is an example pipeline configuration to load using reseqtrack/scripts/pipeline/load_pipeline_from_conf.pl

[alignment]
table_name=sample
config_module=ReseqTrack::Hive::PipeConfig::Alignment_conf
config_options=-root_output_dir /path/to/dir
config_options=-reference /path/to/human.fa
config_options=-known_indels_vcf /nfs/1000g-archive/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz
config_options=-known_snps_vcf /nfs/1000g-archive/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz
config_options=-realign_intervals_file /path/to/realigner.intervals
config_options=-realign_level 1
config_options=-recalibrate_level 1
config_options=-require_experiment_columns library_strategy=WGS
config_options=-sample_attributes POPULATION
config_options=-final_output_layout '#POPULATION#/#sample_alias#/alignment'

  Options that MUST be specified in the pipeline.config_options table/column of your ReseqTrack database:

      -reference, fasta file of your reference genome.  Should be indexed for bwa and should have a .fai and .dict

  Options that have defaults but you will often want to set them in your pipeline.cofig_options table/column:

      -seeding_module, (default is ReseqTrack::Hive::PipeSeed::BasePipeSeed) override this with a project-specific module
      -seeding_options, hashref passed to the seeding module.  Override the defaults only if using a different seeding module.

      -name_file_module, (default is ReseqTrack::Hive::NameFile::BaseNameFile) override this with a project-specific module. Controls how your output bam file is named.
      -name_file_method, (default is 'basic'), controls which subroutine of your name_file_module is used to name your bam file.
      -final_output_dir, (default is your root_output_dir) the root output directory for your final bam file
      -final_output_layout, (default is '#sample_source_id#/alignment') the sub directory for your final bam file
      -name_file_params. a hash ref, passed to your name_file_module.  Change the default values to control how your final output file is named.
            e.g. -name_file_params add_datestamp=0 -new_basename='#sample_alias#.bwa-mem.hg38'

      -RGSM, read group sample appearing in bam file.  Default is '#sample_source_id#'.  Change this to '#sample_alias#' if you want the sample_alias to appear in your bam header
      -RGPU, read group platform unit.  Default is '#sample_run_id#'.  Change this to '#run_alias#' if you want the run_alias to appear in your bam header

      -root_output_dir, (default is your current directory) This is where working files go, i.e. not necessarily the final resting place of your output bam
      -type_fastq, type of fastq files to look for in the reseqtrack database, default FILTERED_FASTQ

      -chunk_max_reads, (default 5000000) controls how fastq files are split up into chunks for parallel alignment

      Recalibration and Realignment options:
      -realign_level, can be 0 (don't realign), 1 (fast realignment around known indels only at lane level), 2 (full realignment at sample level).
      -recalibrate_level, can be 0 (don't recalibrate), 1 (fast recalibration at lane level, e.g. 1000genomes), 2 (slower recalibration at sample level for better accuracy)
      -run_calmd, can be 0 (don't don't run calmd), 1 (do run it, the default)
      -known_indels_vcf, used for indel realignment (default undefined).  Optional if realign_level=2; mandatory if realign_level=1.
      -known_snps_vcf, used for recalibration (default undefined). Mandatory if recalibrate_level != 0.
      -realign_intervals_file, should be given if realign_level=1.  Can be generated using gatk RealignerTargetCreator
      -recalibration_chromosomes, run recalibration_run_level for specified chromosomes; default is for GRCh38 recalibration_chromosomes=>join(',', (map {'chr'.$_} (1..22, 'X', 'Y', 'M'))). 
            For earlier reference genomes, can be recalibration_chromosomes=>join(',', 1..22, 'X', 'Y', 'MT')

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
      -biobambam_dir => (default /nfs/1000g-work/G1K/work/bin/biobambam/bin)

      -sample_columns, default is ['sample_id', 'sample_source_id', 'sample_alias'].
      -run_columns, default is ['run_source_id', 'center_name', 'run_alias'],
      -study_columns, default is ['study_source_id']
      -experiment_columns, default is ['instrument_platform', 'paired_nominal_length'],
      -sample attributes, -run_attributes, -experiment_attributes, study_attributes, default is [] for each one.
            These parameters define what meta information parameters are added to the flow of information around the hive pipeline
            Add to these arrays if your pipeline uses any extra meta information, e.g. when naming the final output files.
            e.g. for 1000genomes project you might want -sample_attributes POPULATION

      -require_experiment_columns, default is { instrument_platform => ['ILLUMINA'] }
      -require_run_columns, default is { status => ['public', 'private'], }
      -require_run_attributes, -require_experiment_attributes, -require_study_attributes, -require_sample_attributes, default is {} for each one.
      -exclude_run_attributes, -exclude_experiment_attributes, -exclude_study_attributes, -exclude_sample_attributes, default is {} for each one.
      -require_study_columns, -require_sample_columns, default is {} for each one.
      -exclude_sample_columns, default is {}
            Use these hashrefs to control what runs are used for alignment
            e.g. -require_experiment_columns library_strategy=WGS
            e.g. -exclude_sample_attributes SEX=male
            e.g. -require_sample_columns tax_id="#expr([9925,9923])expr#'    (i.e. goats only)

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
        %{$self->SUPER::default_options()},

        'pipeline_name'               => 'align',

        seeding_module                => 'ReseqTrack::Hive::PipeSeed::BasePipeSeed',
        seeding_options               => {
            output_columns     => $self->o('sample_columns'),
            output_attributes  => $self->o('sample_attributes'),
            require_columns    => $self->o('require_sample_columns'),
            exclude_columns    => $self->o('exclude_sample_columns'),
            require_attributes => $self->o('require_sample_attributes'),
            exclude_attributes => $self->o('exclude_sample_attributes'),
        },

        'chunk_max_reads'             => 5000000,
        'type_fastq'                  => 'FILTERED_FASTQ',
        'split_exe'                   => $self->o('ENV', 'RESEQTRACK') . '/c_code/split/split',
        'validate_bam_exe'            => $self->o('ENV', 'RESEQTRACK') . '/c_code/validate_bam/validate_bam',
        'bwa_exe'                     => '/nfs/1000g-work/G1K/work/bin/bwa/bwa',
        'samtools_exe'                => '/nfs/1000g-work/G1K/work/bin/samtools/samtools',
        'cramtools_dir'               => '/nfs/1000g-work/G1K/work/bin/crammer',
        'cramtools_jar_file'          => 'cramtools-2.1.jar',
        'squeeze_exe'                 => '/nfs/1000g-work/G1K/work/bin/bamUtil/bin/bam',
        'gatk_dir'                    => '/nfs/1000g-work/G1K/work/bin/gatk/',
        'picard_dir'                  => '/nfs/1000g-work/G1K/work/bin/picard',
        'biobambam_dir'               => '/nfs/1000g-work/G1K/work/bin/biobambam/bin',
        'known_indels_vcf'            => [],
        'known_snps_vcf'              => undef,
        'realign_intervals_file'      => undef,
        'max_paired_length'           => undef,

        #various options for overriding defaults in reheadering bam
        'dict_file'                   => undef,
        'reference_uri'               => undef,
        'ref_assembly'                => undef,
        'ref_species'                 => undef,
        'header_lines_file'           => undef,

        #recalibration_chromosomes => join(',', (map {'chr'.$_} (1..22, 'X', 'Y', 'M'))),
        recalibration_chromosomes     => undef,

        'realign_level'               => 0,
        'recalibrate_level'           => 2,
        'run_calmd'                   => 1,

        'bam_type'                    => undef,
        'bai_type'                    => undef,
        'bas_type'                    => undef,

        'cram_type'                   => undef,
        'crai_type'                   => undef,

        'bwa_algorithm'               => 'mem',
        'bwa_options'                 => { algorithm => $self->o('bwa_algorithm'), threads => 3},

        gatk_threads                  => 3,

        'RGSM'                        => '#sample_source_id#',
        'RGPU'                        => '#run_source_id#',

        'sample_attributes'           => [],
        'sample_columns'              => [ 'sample_id', 'sample_source_id', 'sample_alias' ],
        'run_attributes'              => [],
        'run_columns'                 => [ 'run_source_id', 'center_name', 'run_alias' ],
        'study_attributes'            => [],
        'study_columns'               => [ 'study_source_id' ],
        'experiment_attributes'       => [],
        'experiment_columns'          => [ 'instrument_platform', 'paired_nominal_length' ],

        require_run_attributes        => {},
        require_experiment_attributes => {},
        require_study_attributes      => {},
        require_sample_attributes     => {},
        exclude_run_attributes        => {},
        exclude_experiment_attributes => {},
        exclude_study_attributes      => {},
        exclude_sample_attributes     => {},
        require_experiment_columns    => { instrument_platform => [ 'ILLUMINA' ], },
        require_run_columns           => { status => [ 'public', 'private' ], },
        require_study_columns         => {},
        require_sample_columns        => {},
        exclude_sample_columns        => {},


        final_output_dir              => $self->o('root_output_dir'),
        final_output_layout           => '#sample_source_id#/alignment',
        name_file_module              => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method              => 'basic',
        name_file_params              => {
            new_dir       => '#final_output_dir#/#final_output_layout#',
            new_basename  => '#sample_source_id#.bwa',
            add_datestamp => 1,
            suffix        => '.bam',
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

        dir_label_params => [ "sample_source_id", "library_name", "run_source_id", "chunk_label" ],
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes}, # inherit 'default' from the parent class
        '200Mb' => { 'LSF' => '-C0 -M200 -q ' . $self->o('lsf_queue') . ' -R"select[mem>200] rusage[mem=200] select[' . $self->o('lsf_resource') . ']"' },
        '500Mb' => { 'LSF' => '-C0 -M500 -q ' . $self->o('lsf_queue') . ' -R"select[mem>500] rusage[mem=500] select[' . $self->o('lsf_resource') . ']"' },
        '800Mb' => { 'LSF' => '-C0 -M800 -q ' . $self->o('lsf_queue') . ' -R"select[mem>800] rusage[mem=800] select[' . $self->o('lsf_resource') . ']"' },
        '1Gb'   => { 'LSF' => '-C0 -M1000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>1000] rusage[mem=1000] select[' . $self->o('lsf_resource') . ']"' },
        '2Gb'   => { 'LSF' => '-C0 -M2000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>2000] rusage[mem=2000] select[' . $self->o('lsf_resource') . ']"' },
        '3Gb'   => { 'LSF' => '-C0 -M3000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>3000] rusage[mem=3000] select[' . $self->o('lsf_resource') . ']"' },
        '3Gb3cpus'   => { 'LSF' => '-n 3 -C0 -M3000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>3000] rusage[mem=3000] select[' . $self->o('lsf_resource') . ']"' },
        '4Gb'   => { 'LSF' => '-C0 -M4000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>4000] rusage[mem=4000] select[' . $self->o('lsf_resource') . ']"' },
        '5Gb'   => { 'LSF' => '-C0 -M5000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>5000] rusage[mem=5000] select[' . $self->o('lsf_resource') . ']"' },
        '5Gb3cpus' => { 'LSF' => '-n 3 -C0 -M5000 -q '.$self->o('lsf_queue').' -R"select[mem>5000] rusage[mem=5000] select[' . $self->o('lsf_resource') . ']"' },
        '6Gb'   => { 'LSF' => '-C0 -M6000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>6000] rusage[mem=6000] select[' . $self->o('lsf_resource') . ']"' },
        '8Gb'   => { 'LSF' => '-C0 -M8000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>8000] rusage[mem=8000] select[' . $self->o('lsf_resource') . ']"' },
        '8Gb3cpus' => { 'LSF' => '-n 3 -C0 -M8000 -q '.$self->o('lsf_queue').' -R"select[mem>8000] rusage[mem=8000] select[' . $self->o('lsf_resource') . ']"' },
        '10Gb'  => { 'LSF' => '-C0 -M10000 -q ' . $self->o('lsf_queue') . ' -R"select[mem>10000] rusage[mem=10000] select[' . $self->o('lsf_resource') . ']"' },
    };
}


sub pipeline_analyses {
    my ($self) = @_;

    my @analyses;
    push(@analyses, {
        -logic_name  => 'get_seeds',
        -module      => 'ReseqTrack::Hive::Process::SeedFactory',
        -meadow_type => 'LOCAL',
        -parameters  => {
            seeding_module  => $self->o('seeding_module'),
            seeding_options => $self->o('seeding_options'),
        },
        -flow_into   => {
            2 => [ 'libraries_factory' ],
        },
    });

    push(@analyses, {
        -logic_name  => 'libraries_factory',
        -module      => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
        -meadow_type => 'LOCAL',
        -parameters  => {
            factory_type                  => 'library',
            require_experiment_columns    => $self->o('require_experiment_columns'),
            require_study_columns         => $self->o('require_study_columns'),
            require_experiment_attributes => $self->o('require_experiment_attributes'),
            require_study_attributes      => $self->o('require_study_attributes'),
            exclude_experiment_attributes => $self->o('exclude_experiment_attributes'),
            exclude_study_attributes      => $self->o('exclude_study_attributes'),
        },
        -flow_into   => {
            '2->A' => [ 'runs_factory' ],
            'A->1' => [ 'decide_merge_libraries' ],
        },
    });

    push(@analyses, {
        -logic_name  => 'runs_factory',
        -module      => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
        -meadow_type => 'LOCAL',
        -parameters  => {
            factory_type                  => 'run',
            require_experiment_columns    => $self->o('require_experiment_columns'),
            require_study_columns         => $self->o('require_study_columns'),
            require_run_columns           => $self->o('require_run_columns'),
            require_experiment_attributes => $self->o('require_experiment_attributes'),
            require_study_attributes      => $self->o('require_study_attributes'),
            require_run_attributes        => $self->o('require_run_attributes'),
            exclude_experiment_attributes => $self->o('exclude_experiment_attributes'),
            exclude_study_attributes      => $self->o('exclude_study_attributes'),
            exclude_run_attributes        => $self->o('exclude_run_attributes'),

            output_run_columns            => $self->o('run_columns'),
            output_study_columns          => $self->o('study_columns'),
            output_experiment_columns     => $self->o('experiment_columns'),
            output_run_attributes         => $self->o('run_attributes'),
            output_study_attributes       => $self->o('study_attributes'),
            output_experiment_attributes  => $self->o('experiment_attributes'),
        },
        -flow_into   => {
            '2->A' => [ 'find_source_fastqs' ],
            'A->1' => [ 'decide_mark_duplicates' ],
        },
    });

    push(@analyses, {
        -logic_name  => 'find_source_fastqs',
        -module      => 'ReseqTrack::Hive::Process::ImportCollection',
        -meadow_type => 'LOCAL',
        -parameters  => {
            collection_type    => $self->o('type_fastq'),
            collection_name    => '#run_source_id#',
            output_param       => 'fastq',
            reseqtrack_options => {
                flows_do_count_param => 'fastq',
                flows_do_count       => { 1 => '1+', },
            }
        },
        -flow_into   => {
            1 => [ 'split_fastq', ':////accu?fastq=[]' ],
        },
    });

    push(@analyses, {
        -logic_name        => 'split_fastq',
        -module            => 'ReseqTrack::Hive::Process::SplitFastq',
        -parameters        => {
            program_file => $self->o('split_exe'),
            max_reads    => $self->o('chunk_max_reads'),
        },
        -rc_name           => '200Mb',
        -analysis_capacity => 4,
        -hive_capacity     => 200,
        -flow_into         => {
            '2->A' => { 'bwa' => { 'chunk_label' => '#run_source_id#.#chunk#', 'fastq' => '#fastq#' } },
            'A->1' => [ 'decide_merge_chunks' ],
        }
    });

    push(@analyses, {
        -logic_name    => 'bwa',
        -module        => 'ReseqTrack::Hive::Process::BWA',
        -parameters    => {
            program_file       => $self->o('bwa_exe'),
            samtools           => $self->o('samtools_exe'),
            reference          => $self->o('reference'),
            options            => $self->o('bwa_options'),
            RGSM               => $self->o('RGSM'),
            RGPU               => $self->o('RGPU'),
            reseqtrack_options => {
                delete_param => 'fastq',
            },
        },
        -rc_name       => '8Gb3cpus', # Note the 'hardened' version of BWA may need 8Gb RAM or more
        -hive_capacity => 100,
        -flow_into     => {
            1 => [ 'sort_chunks' ],
        },
    });

=head    
    push(@analyses, {
            -logic_name => 'sort_chunks',
            -module        => 'ReseqTrack::Hive::Process::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                bwa_algorithm => $self->o('bwa_algorithm'),
                command => '#expr(#bwa_algorithm#=="sw" ? "add_or_replace_read_groups" : "fix_mate")expr#',
                options => {read_group_fields => {
                  ID => '#run_source_id#',
                  LB => '#library_name#',
                  PL => '#instrument_platform#',
                  PU => $self->o('RGPU'),
                  SM => $self->o('RGSM'),
                  CN => '#center_name#',
                  DS => '#study_source_id#',
                  PI => '#paired_nominal_length#',
                }, },
                create_index => 1,
                jvm_args => '-Xmx2g',
                reseqtrack_options => {
                  delete_param => 'bam',
                },
            },
            -rc_name => '2Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [ ':////accu?bam=[]', ':////accu?bai=[]']
            },
      });
=cut      

    ### As BioBamBam seems like not having a function called add_or_replace_read_groups, cannot do that for bwa-sw. When bwa-sw is used, one needs to switch to
    ### the picard sort (above) for sorting.
    ### the resource requirement is 2G for BioBamBam Bamsort as well; no saving in memory comparing with Picard

    push(@analyses, {
        -logic_name    => 'sort_chunks',
        -module        => 'ReseqTrack::Hive::Process::RunBioBamBam',
        -parameters    => {
            biobambam_dir      => $self->o('biobambam_dir'),
            command            => "bamsort",
            options            => { fixmates => 1 },
            create_index       => 1,
            reseqtrack_options => {
                delete_param => 'bam',
            },
        },
        -rc_name       => '2Gb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ ':////accu?bam=[]', ':////accu?bai=[]' ]
        },
    });

    push(@analyses, {
        -logic_name  => 'decide_merge_chunks',
        -module      => 'ReseqTrack::Hive::Process::BaseProcess',
        -meadow_type => 'LOCAL',
        -parameters  => {
            realign_level      => $self->o('realign_level'),
            recalibrate_level  => $self->o('recalibrate_level'),
            run_calmd          => $self->o('run_calmd'),
            reseqtrack_options => {
                flows_non_factory    => {
                    1 => '#expr(#realign_level#==1)expr#',
                    2 => '#expr(#realign_level#==1)expr#',
                    3 => '#expr(#recalibrate_level#==1 && (#realign_level#==2 || (#realign_level#==0 && #run_calmd#==0)))expr#',
                    4 => '#expr(#recalibrate_level#==1 && (#realign_level#==2 || (#realign_level#==0 && #run_calmd#==0)))expr#',
                    5 => '#expr(#realign_level#==0 && #run_calmd#==1)expr#',
                    6 => '#expr(#realign_level#==0 && #run_calmd#==1)expr#',
                    7 => '#expr(#realign_level#==0 && #recalibrate_level#==0 && #run_calmd#==0)expr#',
                    8 => '#expr(#realign_level#==0 && #recalibrate_level#==0 && #run_calmd#==0)expr#',
                    9 => '#expr(#realign_level#==2 && #recalibrate_level#!=1)expr#',
                },
                flows_do_count_param => 'bam',
                flows_do_count       => {
                    1 => '1+',
                    2 => '2+',
                    3 => '1+',
                    4 => '2+',
                    5 => '1+',
                    6 => '2+',
                    7 => '1+',
                    8 => '2+',
                    9 => '1+',
                },
            },
        },
        -flow_into   => {
            '2->A' => [ 'merge_chunks' ],
            'A->1' => [ 'realign_knowns_only' ],
            '4->B' => [ 'merge_chunks' ],
            'B->3' => [ 'recalibrate_run_level' ],
            '6->C' => [ 'merge_chunks' ],
            'C->5' => [ 'calmd_run_level' ],
            '8->D' => [ 'merge_chunks' ],
            'D->7' => [ 'tag_strip_run_level' ],
            9      => [ ':////accu?bam=[]', ':////accu?bai=[]' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'merge_chunks',
        -module        => 'ReseqTrack::Hive::Process::RunBioBamBam',
        -parameters    => {
            biobambam_dir      => $self->o('biobambam_dir'),
            command            => 'bammerge',
            create_index       => 1,
            reseqtrack_options => {
                delete_param => [ 'bam', 'bai' ],
            }
        },
        -rc_name       => '500Mb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ ':////accu?bam=[]', ':////accu?bai=[]' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'realign_knowns_only',
        -module        => 'ReseqTrack::Hive::Process::RunBamImprovement',
        -parameters    => {
            command             => 'realign',
            reference           => $self->o('reference'),
            gatk_dir            => $self->o('gatk_dir'),
            known_sites_vcf     => $self->o('known_indels_vcf'),
            intervals_file      => $self->o('realign_intervals_file'),
            gatk_module_options => {
                knowns_only => 1,
                threads     => $self->o('gatk_threads') },
            recalibrate_level   => $self->o('recalibrate_level'),
            java_exe            => $self->o('java_exe'),
            run_calmd           => $self->o('run_calmd'),
            reseqtrack_options  => {
                delete_param      => [ 'bam', 'bai' ],
                flows_non_factory => {
                    1 => '#expr(#run_calmd#==1)expr#',
                    2 => '#expr(#run_calmd#==0 && #recalibrate_level#==1)expr#',
                    3 => '#expr(#run_calmd#==0 && #recalibrate_level#==0)expr#',
                    4 => '#expr(#run_calmd#==0 && #recalibrate_level#==2)expr#',
                }
            }
        },
        -rc_name       => '5Gb3cpus',
        -hive_capacity => 100,
        -flow_into     => {
            1 => [ 'calmd_run_level' ],
            2 => [ 'index_recalibrate_run_level' ],
            3 => [ 'tag_strip_run_level' ],
            4 => [ ':////accu?bam=[]', ':////accu?bai=[]' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'calmd_run_level',
        -module        => 'ReseqTrack::Hive::Process::RunSamtools',
        -parameters    => {
            program_file       => $self->o('samtools_exe'),
            command            => 'calmd',
            reference          => $self->o('reference'),
            samtools_options   => { input_sort_status => 'c' },
            recalibrate_level  => $self->o('recalibrate_level'),
            reseqtrack_options => {
                delete_param      => [ 'bam' ],
                flows_non_factory => {
                    1 => '#expr(#recalibrate_level#==1)expr#',
                    2 => '#expr(#recalibrate_level#==0)expr#',
                    3 => '#expr(#recalibrate_level#==2)expr#',
                }
            }
        },
        -rc_name       => '2Gb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ 'index_recalibrate_run_level' ],
            2 => [ 'tag_strip_run_level' ],
            3 => [ ':////accu?bam=[]', ':////accu?bai=[]' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'index_recalibrate_run_level',
        -module        => 'ReseqTrack::Hive::Process::RunSamtools',
        -parameters    => {
            program_file => $self->o('samtools_exe'),
            command      => 'index',
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ 'recalibrate_run_level' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'recalibrate_run_level',
        -module        => 'ReseqTrack::Hive::Process::RunBamImprovement',
        -parameters    => {
            command             => 'recalibrate',
            reference           => $self->o('reference'),
            gatk_dir            => $self->o('gatk_dir'),
            jvm_args            => '-Xmx2g',
            java_exe            => $self->o('java_exe'),
            known_sites_vcf     => $self->o('known_snps_vcf'),
            realign_level       => $self->o('realign_level'),
            gatk_module_options => {
                intervals   => $self->o('recalibration_chromosomes'),
                # This tool requires that parallelism is set by using the -nct option.
                threads     => 1 },
            reseqtrack_options  => {
                delete_param      => [ 'bam', 'bai' ],
                flows_non_factory => {
                    1 => '#expr(#realign_level#!=2)expr#',
                    2 => '#expr(#realign_level#==2)expr#',
                }
            }
        },
        -rc_name       => '3Gb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ 'tag_strip_run_level' ],
            2 => [ ':////accu?bam=[]', ':////accu?bai=[]' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'tag_strip_run_level',
        -module        => 'ReseqTrack::Hive::Process::RunSqueezeBam',
        -parameters    => {
            program_file       => $self->o('squeeze_exe'),
            'rm_OQ_fields'     => 1,
            'rm_tag_types'     => [ 'XM:i', 'XG:i', 'XO:i', 'pa:f', 'BD:Z', 'BI:Z' ],
            reseqtrack_options => {
                delete_param => [ 'bam' ],
            }
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ ':////accu?bam=[]', ':////accu?bai=[]' ],
        },
    });

    ## remove_tag_type pa:f is most likely from biobambam, it messes up with bam validate, the step that generates bas
    push(@analyses, {
        -logic_name  => 'decide_mark_duplicates',
        -module      => 'ReseqTrack::Hive::Process::BaseProcess',
        -meadow_type => 'LOCAL',
        -parameters  => {
            reseqtrack_options => {
                denestify            => 'bam',
                flows_do_count_param => 'bam',
                flows_do_count       => { 1 => '1+', },
            }
        },
        -flow_into   => {
            1 => [ 'mark_duplicates', ':////accu?fastq=[]' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'mark_duplicates',
        -module        => 'ReseqTrack::Hive::Process::RunBioBamBam',
        -parameters    => {
            biobambam_dir      => $self->o('biobambam_dir'),
            command            => 'bammarkduplicates',
            create_index       => 1,
            reseqtrack_options => {
                delete_param => [ 'bam', 'bai' ],
                denestify    => [ 'bam', 'bai' ],
            }
        },
        -rc_name       => '800Mb',
        -hive_capacity => 100,
        -flow_into     => {
            1 => [ ':////accu?bam=[]', ':////accu?bai=[]' ],
        },
    });

    push(@analyses, {
        -logic_name  => 'decide_merge_libraries',
        -module      => 'ReseqTrack::Hive::Process::BaseProcess',
        -meadow_type => 'LOCAL',
        -parameters  => {
            realign_level      => $self->o('realign_level'),
            recalibrate_level  => $self->o('recalibrate_level'),
            reseqtrack_options => {
                flows_non_factory    => {
                    1 => '#expr(#realign_level#==2)expr#',
                    2 => '#expr(#realign_level#==2)expr#',
                    3 => '#expr(#recalibrate_level#==2 && #realign_level#!=2)expr#',
                    4 => '#expr(#recalibrate_level#==2 && #realign_level#!=2)expr#',
                    5 => '#expr(#recalibrate_level#!=2 && #realign_level#!=2)expr#',
                    6 => '#expr(#recalibrate_level#!=2 && #realign_level#!=2)expr#',
                    7 => 1,
                },
                flows_do_count_param => 'bam',
                flows_do_count       => {
                    1 => '1+',
                    2 => '2+',
                    3 => '1+',
                    4 => '2+',
                    5 => '1+',
                    6 => '2+',
                    7 => '0',
                },
            },
        },
        -flow_into   => {
            '2->A' => [ 'merge_libraries' ],
            'A->1' => [ 'realign_full' ],
            '4->B' => [ 'merge_libraries' ],
            'B->3' => [ 'recalibrate_sample_level' ],
            '6->C' => [ 'merge_libraries' ],
            'C->5' => [ 'reheader' ],
            '7'    => [ 'mark_seed_futile' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'merge_libraries',
        -module        => 'ReseqTrack::Hive::Process::RunBioBamBam',
        -parameters    => {
            biobambam_dir      => $self->o('biobambam_dir'),
            command            => 'bammerge',
            create_index       => 1,
            reseqtrack_options => {
                delete_param => [ 'bam', 'bai' ],
            },
        },
        -rc_name       => '500Mb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ ':////accu?bam=[]', ':////accu?bai=[]' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'realign_full',
        -module        => 'ReseqTrack::Hive::Process::RunBamImprovement',
        -parameters    => {
            command             => 'realign',
            reference           => $self->o('reference'),
            gatk_dir            => $self->o('gatk_dir'),
            java_exe            => $self->o('java_exe'),
            known_sites_vcf     => $self->o('known_indels_vcf'),
            gatk_module_options => {
                knowns_only => 0,
                threads => $self->o('gatk_threads') },
            recalibrate_level   => $self->o('recalibrate_level'),
            run_calmd           => $self->o('run_calmd'),
            reseqtrack_options  => {
                delete_param      => [ 'bam', 'bai' ],
                flows_non_factory => {
                    1 => '#expr(#run_calmd#==1)expr#',
                    2 => '#expr(#run_calmd#==0 && #recalibrate_level#==2)expr#',
                    3 => '#expr(#run_calmd#==0 && #recalibrate_level#!=2)expr#',
                }
            },
        },
        -rc_name       => '5Gb3cpus',
        -hive_capacity => 100,
        -flow_into     => {
            1 => [ 'calmd_sample_level' ],
            2 => [ 'index_recalibrate_sample_level' ],
            3 => [ 'tag_strip_sample_level' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'calmd_sample_level',
        -module        => 'ReseqTrack::Hive::Process::RunSamtools',
        -parameters    => {
            program_file       => $self->o('samtools_exe'),
            command            => 'calmd',
            reference          => $self->o('reference'),
            samtools_options   => { input_sort_status => 'c' },
            recalibrate_level  => $self->o('recalibrate_level'),
            reseqtrack_options => {
                delete_param      => [ 'bam' ],
                flows_non_factory => {
                    1 => '#expr(#recalibrate_level#==2)expr#',
                    2 => '#expr(#recalibrate_level#!=2)expr#',
                }
            }
        },
        -rc_name       => '2Gb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ 'index_recalibrate_sample_level' ],
            2 => [ 'tag_strip_sample_level' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'index_recalibrate_sample_level',
        -module        => 'ReseqTrack::Hive::Process::RunSamtools',
        -parameters    => {
            program_file => $self->o('samtools_exe'),
            command      => 'index',
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ 'recalibrate_sample_level' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'recalibrate_sample_level',
        -module        => 'ReseqTrack::Hive::Process::RunBamImprovement',
        -parameters    => {
            command            => 'recalibrate',
            reference          => $self->o('reference'),
            gatk_dir           => $self->o('gatk_dir'),
            gatk_module_options => { threads => $self->o('gatk_threads') },
            jvm_args           => '-Xmx2g',
            java_exe           => $self->o('java_exe'),
            known_sites_vcf    => $self->o('known_snps_vcf'),
            reseqtrack_options => {
                delete_param => [ 'bam', 'bai' ],
            },
        },
        -rc_name       => '3Gb3cpus',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ 'tag_strip_sample_level' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'tag_strip_sample_level',
        -module        => 'ReseqTrack::Hive::Process::RunSqueezeBam',
        -parameters    => {
            program_file       => $self->o('squeeze_exe'),
            'rm_OQ_fields'     => 1,
            'rm_tag_types'     => [ 'XM:i', 'XG:i', 'XO:i', 'pa:f', 'BD:Z', 'BI:Z' ],
            reseqtrack_options => {
                delete_param => [ 'bam' ],
            },
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ 'reheader' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'reheader',
        -module        => 'ReseqTrack::Hive::Process::ReheaderBam',
        -parameters    => {
            'samtools'          => $self->o('samtools_exe'),
            'header_lines_file' => $self->o('header_lines_file'),
            'dict_file'         => $self->o('dict_file'),
            'reference'         => $self->o('reference'),
            'SQ_assembly'       => $self->o('ref_assembly'),
            'SQ_species'        => $self->o('ref_species'),
            'SQ_uri'            => $self->o('reference_uri'),
            reseqtrack_options  => {
                denestify    => 'fastq',
                delete_param => [ 'bam', 'bai' ],
            },
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => [ 'final_index' ],
        },
    });

    push(@analyses, {
        -logic_name    => 'final_index',
        -module        => 'ReseqTrack::Hive::Process::RunSamtools',
        -parameters    => {
            program_file => $self->o('samtools_exe'),
            command      => 'index',
        },
        -flow_into     => { 1 => [ 'store_bam' ] },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
    });

    push(@analyses, {
        -logic_name    => 'store_bam',
        -module        => 'ReseqTrack::Hive::Process::LoadFile',
        -parameters    => {
            type                => $self->o('bam_type'),
            file                => '#bam#',
            name_file_module    => $self->o('name_file_module'),
            name_file_method    => $self->o('name_file_method'),
            name_file_params    => $self->o('name_file_params'),
            final_output_dir    => $self->o('final_output_dir'),
            final_output_layout => $self->o('final_output_layout'),
            clobber             => 1,
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => { 1 => { 'store_bai' => { 'bam' => '#file#' } } },
    });

    push(@analyses, {
        -logic_name    => 'store_bai',
        -module        => 'ReseqTrack::Hive::Process::LoadFile',
        -parameters    => {
            type             => $self->o('bai_type'),
            file             => '#bai#',
            name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
            name_file_method => 'basic',
            name_file_params => { new_full_path => '#bam#.bai' },
            clobber          => 1,
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => { 1 => { 'validate' => { 'bai' => '#file#' } } },
    });

    push(@analyses, {
        -logic_name    => 'validate',
        -module        => 'ReseqTrack::Hive::Process::RunValidateBam',
        -parameters    => {
            'program_file' => $self->o('validate_bam_exe'),
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => { 1 => [ 'store_bas' ] },
    });

    push(@analyses, {
        -logic_name  => 'store_bas',
        -module      => 'ReseqTrack::Hive::Process::LoadFile',
        -meadow_type => 'LOCAL',
        -parameters  => {
            type             => $self->o('bas_type'),
            file             => '#bas#',
            name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
            name_file_method => 'basic',
            name_file_params => { new_full_path => '#bam#.bas' },
            clobber          => 1,
        },
        -flow_into   => { 1 => [ 'bam_to_cram' ] },
    });

    push(@analyses, {
        -logic_name    => 'bam_to_cram',
        -module        => 'ReseqTrack::Hive::Process::RunCramtools',
        -parameters    => {
            bam                => '#bam#',
            java_exe           => $self->o('java_exe'),
            jvm_args           => '-Xmx4g',
            program_file       => $self->o('cramtools_dir'),
            cramtools_jar_file => $self->o('cramtools_jar_file'),
            reference          => $self->o('reference'),
            cramtools_options  => { 'preserve-read-names' => 1, 'capture-all-tags' => 1, 'ignore-tags' => 'OQ:CQ:BQ', 'lossy-quality-score-spec' => '\'*8\'' },
        },
        -rc_name       => '4Gb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => { 'store_cram' => { 'file' => '#cram#' } },
        },
    });

    push(@analyses, {
        -logic_name    => 'store_cram',
        -module        => 'ReseqTrack::Hive::Process::LoadFile',
        -parameters    => {
            type                => $self->o('cram_type'),
            file                => '#cram#',
            name_file_module    => $self->o('name_file_module'),
            name_file_method    => $self->o('name_file_method'),
            name_file_params    => { new_full_path => '#bam#.cram' },
            final_output_dir    => $self->o('final_output_dir'),
            final_output_layout => $self->o('final_output_layout'),
            clobber             => 1,
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => { 'index_cram' => { 'bam' => '#file#' } },
        },
    });

    push(@analyses, {
        -logic_name    => 'index_cram',
        -module        => 'ReseqTrack::Hive::Process::RunCramtools',
        -parameters    => {
            bam                => '#file#',
            java_exe           => $self->o('java_exe'),
            jvm_args           => '-Xmx2g',
            program_file       => $self->o('cramtools_dir'),
            cramtools_jar_file => $self->o('cramtools_jar_file'),
        },
        -rc_name       => '2Gb',
        -hive_capacity => 200,
        -flow_into     => {
            1 => { 'store_crai' => { 'file' => '#crai#' } },
        },
    });

    push(@analyses, {
        -logic_name    => 'store_crai',
        -module        => 'ReseqTrack::Hive::Process::LoadFile',
        -parameters    => {
            type             => $self->o('crai_type'),
            name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
            name_file_method => 'basic',
            name_file_params => { new_full_path => '#file#' },
            clobber          => 1,
        },
        -rc_name       => '200Mb',
        -hive_capacity => 200,
        -flow_into     => { 1 => { 'mark_seed_complete' => { 'crai' => '#file#' } } },
    });

    push(@analyses, {
        -logic_name  => 'mark_seed_complete',
        -module      => 'ReseqTrack::Hive::Process::UpdateSeed',
        -parameters  => {
            is_complete => 1,
        },
        -meadow_type => 'LOCAL',
    });

    push(@analyses, {
        -logic_name  => 'mark_seed_futile',
        -module      => 'ReseqTrack::Hive::Process::UpdateSeed',
        -parameters  => {
            is_futile => 1,
        },
        -meadow_type => 'LOCAL',
    });
    return \@analyses;
}

1;

