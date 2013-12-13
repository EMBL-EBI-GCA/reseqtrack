=head1 NAME

 ReseqTrack::Hive::PipeConfig::VQSR_conf

=head1 SYNOPSIS

  Pipeline must be seeded by the collection table of a ReseqTrack database (collection of bam files)
  Vcf files will be created for each collection (possibly for SNPs and for indels), and stored in the ReseqTrack database

  Here is an example pipeline configuration to load using reseqtrack/scripts/pipeline/load_pipeline_from_conf.pl

[vqsr]
table_name=collection
config_module=ReseqTrack::Hive::PipeConfig::VQSR_conf
config_options=-root_output_dir /path/to/dir
config_options=-reference /path/to/human.fa
config_options=-snp_alleles /nfs/1000g-work/G1K/work/REFERENCE/snps/grc37_snps/OMNI25_1856/Omni25_genotypes_1856_PASS.vcf.gz
config_options=-indel_alleles /nfs/1000g-archive/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz
config_options=-snp_resources omni='known=true,training=true,truth=true,prior=12.0 /nfs/1000g-work/G1K/work/REFERENCE/snps/grc37_snps/OMNI25_1856/Omni25_genotypes_1856_PASS.vcf.gz'
config_options=-indel_resources mills='known=true,training=true,truth=true,prior=12.0 /nfs/1000g-archive/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz'
config_options=-callgroup_type CALLGROUP_BAM
  
  Options that MUST be specified in the pipeline.config_options table/column of your ReseqTrack database:

      -callgroup_type, type of collection (of bams) in the reseqtrack database used for seeding the pipeline

      -reference, fasta file of your reference genome.  Should have a .fai

      -snp_alleles, path to a vcf file containing sites for recalibration (should be tabixed)
      -indel_alleles, path to a vcf file containing sites for recalibration (should be tabixed)

      -snp_resources, key/val pairs e.g. omni='known=true,training=true,truth=true,prior=12.0 /path/to/file.vcf
      -indel_resources, key/val pairs e.g.  mills='known=true,training=true,truth=true,prior=12.0 /path/to/file.vcf
          Resource files should be tabixed

  Options that have defaults but you will often want to set them in your pipeline.cofig_options table/column:

      -seeding_module, (default is ReseqTrack::Hive::PipeSeed::BasePipeSeed) override this with a project-specific module
      -seeding_options, hashref passed to the seeding module.  Override the defaults only if using a different seeding module.

      -require_collection_columns, -exclude_collection_columns, -require_collection_attributes, -exclude_collection_attributes'
            These are hashrefs, add to these to control which collections are used to seed the pipeline
            e.g. -exclude_collection_columns name=GBR

      -name_file_module, (default is ReseqTrack::Hive::NameFile::BaseNameFile) override this with a project-specific module. Controls how your output bam file is named.
      -name_file_method, (default is 'basic'), controls which subroutine of your name_file_module is used to name your bam file.
      -final_output_dir, (default is your root_output_dir) the root output directory for your final bam file
      -name_file_params. a hash ref, passed to your name_file_module.  Change the default values to control how your final output file is named.

      -root_output_dir, (default is your current directory) This is where working files go, i.e. not necessarily the final resting place of your output vcf
      -target_bed_file, for if you want to do exome calling (default undefined)
      -transpose_window_size, (default 50000000) Controls the size of the region convered by a single transposed bam
      -call_window_size, (default 50000) Controls the size of the region for each individual variant calling job

      Caller options
      -vqsr_snps, boolean, default 1, turn on/off VQSR for SNPs
      -vqsr_indels, boolean, default 1, turn on/off VQSR for indels
      -call_by_gatk_snps_options, key/val pairs, passed on to the CallByGATK module for SNPs
      -call_by_gatk_indels_options, key/val pairs, passed on to the CallByGATK module for indels
      -variant_recalibrator_snps_options, key/val pairs, passed on to the VariantRecalibrator module for SNPs
      -variant_recalibrator_indels_options, key/val pairs, passed on to the VariantRecalibrator module for indels
      -apply_recalibration_snps_options, key/val pairs, passed on to the ApplyRecalibration module for SNPs
      -apply_recalibration_indels_options, key/val pairs, passed on to the ApplyRecalibration module for indels
      -snp_annotations, VQSR modelling is done using these annotations, can be specified multiple times (default is the defaults in the VariantRecalibrator module)
      -indel_annotations, VQSR modelling is done using these annotations, can be specified multiple times (default is the defaults in the VariantRecalibrator module)

      Paths of executables:
      -transpose_exe, (default is to work it out from your environment variable $RESEQTRACK)
      -bgzip_exe, (default /nfs/1000g-work/G1K/work/bin/tabix/bgzip)
      -gatk_dir, (default /nfs/1000g-work/G1K/work/bin/gatk/dist/)

  Options that are required, but will be passed in by reseqtrack/scripts/init_pipeline.pl:

      -pipeline_db -host=???
      -pipeline_db -port=???
      -pipeline_db -user=???
      -dipeline_db -dbname=???
      -reseqtrack_db -host=???
      -reseqtrack_db -user=???
      -reseqtrack_db -port=???
      -reseqtrack_db -pass=???

=cut


package ReseqTrack::Hive::PipeConfig::VQSR_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');


sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },

        'pipeline_name' => 'vc',

        seeding_module => 'ReseqTrack::Hive::PipeSeed::BasePipeSeed',
        seeding_options => {
            output_columns => ['name', 'collection_id'],
            require_columns => $self->o('require_collection_columns'),
            exclude_columns => $self->o('exclude_collection_columns'),
            require_attributes => $self->o('require_collection_attributes'),
            exclude_attributes => $self->o('exclude_collection_attributes'),
          },

        'transpose_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/transpose_bam/transpose_bam',
        'bgzip_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
        'gatk_dir' => '/nfs/1000g-work/G1K/work/bin/gatk/dist/',

        call_by_gatk_snps_options => {
                glm => 'SNP',
                gt_mode => "GENOTYPE_GIVEN_ALLELES",
                out_mode => "EMIT_ALL_SITES",
                alleles => $self->o('snp_alleles'),
              },

        call_by_gatk_indels_options => {
                glm => 'INDEL',
                gt_mode => "GENOTYPE_GIVEN_ALLELES",
                out_mode => "EMIT_ALL_SITES",
                alleles => $self->o('indel_alleles'),
              },

        variant_recalibrator_snps_options => {mode => 'SNP'},
        variant_recalibrator_indels_options => {mode => 'INDEL'},
        apply_recalibration_snps_options => {mode => 'SNP'},
        apply_recalibration_indels_options => {mode => 'INDEL'},

        snp_resources => {},
        indel_resources => {},

        snp_annotations => [], # use defaults
        indel_annotations => [], # use defaults

        'fai' => $self->o('reference') . '.fai',
        'target_bed_file' => undef,

        'call_window_size' => 50000,
        'transpose_window_size' => 10000000,

        'vqsr_snps' => 1,
        'vqsr_indels' => 1,

        'vcf_type' => undef,

        'require_collection_columns' => {'type' => $self->o('callgroup_type')},
        'exclude_collection_columns' => {},
        'require_collection_attributes' => {},
        'exclude_collection_attributes' => {},

        final_output_dir => $self->o('root_output_dir'),
        name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method => 'basic',
        name_file_params => {
            new_dir => $self->o('final_output_dir'),
            new_basename => '#callgroup#.#var_type#.vqsr',
            add_datestamp => 1,
            suffix => ['.vcf.gz'],
          },

    };
}

sub pipeline_create_commands {
    my ($self) = @_;

    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation
    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},          # here we inherit anything from the base class

        fai => $self->o('fai'),

        dir_label_params => ['callgroup', 'region1', 'region2'],
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
                2 => [ 'block_seed_complete' ],
            },
      });
    push(@analyses, {
          -logic_name => 'block_seed_complete',
          -module        => 'ReseqTrack::Hive::Process::BaseProcess',
          -meadow_type=> 'LOCAL',
          -parameters => {
            reseqtrack_options => {
              flows_non_factory => [1,2],
            },
          },
            -flow_into => {
                '2->A' => { 'find_source_bams' => {'callgroup' => '#name#'}},
                'A->1' => [ 'mark_seed_complete' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'find_source_bams',
            -module        => 'ReseqTrack::Hive::Process::ImportCollection',
            -meadow_type => 'LOCAL',
            -parameters    => {
                output_param => 'bam',
                reseqtrack_options => {
                  flows_file_count_param => 'bam',
                  flows_file_count => { 1 => '1+', },
                },
            },
            -flow_into => {
                1 => [ 'regions_factory_1' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'regions_factory_1',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                num_bases => $self->o('transpose_window_size'),
                max_sequences => 2000,
                bed => $self->o('target_bed_file'),
            },
            -flow_into => {
                '2->A' => { 'transpose_bam' => {'region1' => '#callgroup#.#SQ_start#.#bp_start#.#SQ_end#.#bp_end#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                }},
                'A->1' => [ 'decide_mergers'],
            },
      });
    push(@analyses, {
            -logic_name    => 'transpose_bam',
            -module        => 'ReseqTrack::Hive::Process::TransposeBam',
            -parameters    => {
                program_file => $self->o('transpose_exe'),
                create_index => 1,
                bed => $self->o('target_bed_file'),
                region_overlap => 100,
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -flow_into => {
                '1->A' => [ 'regions_factory_2' ],
                'A->1' => [ 'collect_vcf' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'regions_factory_2',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -parameters    => {
                num_bases => $self->o('call_window_size'),
                max_sequences => 1,
                bed => $self->o('target_bed_file'),
                reseqtrack_options => {
                  flows_factory => {
                      1 => $self->o('vqsr_snps'),
                      2 => $self->o('vqsr_indels'),
                      3 => 1,
                  },
                },
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -flow_into => {
                '1' => { 'call_by_gatk_snps' => {'region2' => '#callgroup#.#SQ_start#.#bp_start#.#SQ_end#.#bp_end#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                }},
                '2' => { 'call_by_gatk_indels' => {'region2' => '#callgroup#.#SQ_start#.#bp_start#.#SQ_end#.#bp_end#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                }},
                '3' => [ ':////accu?bp_start=[fan_index]', ':////accu?bp_end=[fan_index]'],
            },
      });
    push(@analyses, {
            -logic_name    => 'collect_vcf',
            -module        => 'ReseqTrack::Hive::Process::BaseProcess',
            -meadow_type => 'LOCAL',
            -parameters    => {
              reseqtrack_options => {
                flows_non_factory => {
                  1 => $self->o('vqsr_snps'),
                  2 => $self->o('vqsr_indels'),
                  3 => 1,
                },
                delete_param => ['bam','bai'],
              },
            },
            -flow_into => {
                1 => [ ':////accu?gatk_snps_vcf=[fan_index]'],
                2 => [ ':////accu?gatk_indels_vcf=[fan_index]'],
                3 => [ ':////accu?bp_start=[fan_index]', ':////accu?bp_end=[fan_index]', ],
            },
      });
    push(@analyses, {
          -logic_name    => 'call_by_gatk_snps',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByGATK',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('call_by_gatk_snps_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => { ':////accu?gatk_snps_vcf=[fan_index]' => {'gatk_snps_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
              -1 => [ 'call_by_gatk_snps_himem' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_gatk_snps_himem',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByGATK',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('call_by_gatk_snps_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => { ':////accu?gatk_snps_vcf=[fan_index]' => {'gatk_snps_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_gatk_indels',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByGATK',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('call_by_gatk_indels_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => { ':////accu?gatk_indels_vcf=[fan_index]' => {'gatk_indels_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
              -1 => [ 'call_by_gatk_indels_himem' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_gatk_indels_himem',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByGATK',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              alleles_file => $self->o('indel_alleles'),
              options => $self->o('call_by_gatk_indels_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => { ':////accu?gatk_indels_vcf=[fan_index]' => {'gatk_indels_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
          },
      });
    push(@analyses, {
          -logic_name => 'decide_mergers',
          -module        => 'ReseqTrack::Hive::Process::BaseProcess',
          -meadow_type=> 'LOCAL',
          -parameters => {
            reseqtrack_options => {
              flows_non_factory => {
                  1 => $self->o('vqsr_snps'),
                  2 => $self->o('vqsr_indels'),
              },
            },
          },
            -flow_into => {
                '1' => { 'merge_vcf_snps' => {'var_type' => 'snps'}},
                '2' => { 'merge_vcf_indels' => {'var_type' => 'indels'}},
            },
      });


    push(@analyses, {
          -logic_name    => 'merge_vcf_snps',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              bgzip => $self->o('bgzip_exe'),
              vcf => '#gatk_snps_vcf#',
              run_tabix => 1,
              reseqtrack_options => {
                denestify => ['vcf','bp_start','bp_end'],
                decode_file_id => 'vcf',
                delete_param => 'vcf',
              },
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
          -flow_into => {
              '1' => [ 'variant_recalibrator_snps' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'merge_vcf_indels',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              bgzip => $self->o('bgzip_exe'),
              vcf => '#gatk_indels_vcf#',
              run_tabix => 1,
              reseqtrack_options => {
                denestify => ['vcf','bp_start','bp_end'],
                decode_file_id => 'vcf',
                delete_param => 'vcf',
              },
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
          -flow_into => {
              '1' => [ 'variant_recalibrator_indels' ],
          },
      });

    push(@analyses, {
          -logic_name    => 'variant_recalibrator_snps',
          -module        => 'ReseqTrack::Hive::Process::RunVariantRecalibrator',
          -parameters    => {
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('variant_recalibrator_snps_options'),
              resources => $self->o('snp_resources'),
              annotations => $self->o('snp_annotations'),
              jvm_args => '-Xmx16g', 
              mode => 'SNP',
          },
          -rc_name => '8Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ 'apply_recalibration_snps' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'variant_recalibrator_indels',
          -module        => 'ReseqTrack::Hive::Process::RunVariantRecalibrator',
          -parameters    => {
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('variant_recalibrator_indels_options'),
              resources => $self->o('indel_resources'),
              annotations => $self->o('indel_annotations'),
              jvm_args => '-Xmx8g', 
              mode => 'INDEL',
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ 'apply_recalibration_indels' ],
          },
      });

    push(@analyses, {
          -logic_name    => 'apply_recalibration_snps',
          -module        => 'ReseqTrack::Hive::Process::RunApplyRecalibration',
          -parameters    => {
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('apply_recalibration_snps_options'),
              reseqtrack_options => {
                delete_param => ['vcf', 'recal_file', 'tbi'],
              },
          },
          -flow_into => { '1' => [ 'store_vcf' ] },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
      });
    push(@analyses, {
          -logic_name    => 'apply_recalibration_indels',
          -module        => 'ReseqTrack::Hive::Process::RunApplyRecalibration',
          -parameters    => {
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('apply_recalibration_indels_options'),
              reseqtrack_options => {
                delete_param => ['vcf', 'recal_file', 'tbi'],
              },
          },
          -flow_into => { '1' => [ 'store_vcf' ] },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
      });
    push(@analyses, {
            -logic_name    => 'store_vcf',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -parameters    => {
              type => $self->o('vcf_type'),
              file => '#vcf#',
              name_file_module => $self->o('name_file_module'),
              name_file_method => $self->o('name_file_method'),
              name_file_params => $self->o('name_file_params'),
            },
            -rc_name => '200Mb',
            -hive_capacity  =>  200,
      });
    push(@analyses, {
            -logic_name    => 'mark_seed_complete',
            -module        => 'ReseqTrack::Hive::Process::UpdateSeed',
            -parameters    => {
              is_complete  => 1,
            },
            -meadow_type => 'LOCAL',
      });
    return \@analyses;
}

1;

