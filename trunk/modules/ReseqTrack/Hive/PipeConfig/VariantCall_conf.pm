=head1 NAME

 ReseqTrack::Hive::PipeConfig::VariantCall_conf

=head1 SYNOPSIS

  Pipeline must be seeded by the collection table of a ReseqTrack database (collection of bam files)
  Multiple vcf files will be created for each collection, and stored in the ReseqTrack database

  Here is an example pipeline configuration to load using reseqtrack/scripts/pipeline/load_pipeline_from_conf.pl

[variant call]
table_name=collection
config_module=ReseqTrack::Hive::PipeConfig::VariantCall_conf
config_options=-root_output_dir /path/to/dir
config_options=-reference /path/to/human.fa
config_options=-coverage_per_sample 12
config_options=-callgroup_type CALLGROUP_BAM
  
  Options that MUST be specified in the pipeline.config_options table/column of your ReseqTrack database:

      -callgroup_type, type of collection (of bams) in the reseqtrack database used for seeding the pipeline
      -reference, fasta file of your reference genome.  Should have a .fai
      -coverage_per_sample, average depth of coverage in each bam file. Passed to CallBySamtools.

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
      -transpose_window_size, (default 10000000) Controls the size of the region convered by a single transposed bam
      -call_window_size, (default 50000) Controls the size of the region for each individual variant calling job

      Caller options
      -call_by_samtools, boolean, default 1, turn on/off calling by samtools mpileup
      -call_by_gatk, boolean, default 1, turn on/off calling by gatk UnifiedGenotyper
      -call_by_freebayes, boolean, default 1, turn on/off calling by freebayes
      -call_by_gatk_options, key/val pairs, passed on to the CallByGATK module
      -call_by_samtools_options, key/val pairs, passed on to the CallBySamtools module
      -call_by_freebayes_options, key/val pairs, passed on to the CallByFreebayes module

      Paths of executables:
      -transpose_exe, (default is to work it out from your environment variable $RESEQTRACK)
      -samtools_exe, (default /nfs/1000g-work/G1K/work/bin/samtools/samtools)
      -bcftools_exe, (default /nfs/1000g-work/G1K/work/bin/samtools/bcftools/bcftools)
      -vcfutils_exe, (default /nfs/1000g-work/G1K/work/bin/samtools/bcftools/vcfutils.pl)
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



package ReseqTrack::Hive::PipeConfig::VariantCall_conf;

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
            output_columns => ['name', 'sample_id'],
            require_columns => $self->o('require_collection_columns'),
            exclude_columns => $self->o('exclude_collection_columns'),
            require_attributes => $self->o('require_collection_attributes'),
            exclude_attributes => $self->o('exclude_collection_attributes'),
          },

        'transpose_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/transpose_bam/transpose_bam',
        'samtools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/samtools',
        'bcftools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/bcftools/bcftools',
        'vcfutils_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/bcftools/vcfutils.pl',
        'bgzip_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
        'gatk_dir' => '/nfs/1000g-work/G1K/work/bin/gatk/dist/',
        'freebayes_exe' => '/nfs/1000g-work/G1K/work/bin/freebayes/bin/freebayes',

        call_by_gatk_options => {}, # use module defaults
        call_by_freebayes_options => {}, # use module defaults

        call_by_samtools_options => {
          depth_of_coverage => '#expr(#num_samples# * #coverage_per_sample#)expr#',
          },

        'fai' => $self->o('reference') . '.fai',
        'target_bed_file' => undef,

        'call_window_size' => 50000,
        'transpose_window_size' => 10000000,

        'call_by_samtools' => 1,
        'call_by_gatk' => 1,
        'call_by_freebayes' => 1,

        'vcf_type' => undef,

        'callgroup_type' => $self->o('callgroup_type'),
        'require_collection_columns' => {'type' => $self->o('callgroup_type')},
        'exclude_collection_columns' => {},
        'require_collection_attributes' => {},
        'exclude_collection_attributes' => {},

        final_output_dir => $self->o('root_output_dir'),
        name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method => 'basic',
        name_file_params => {
            new_dir => $self->o('final_output_dir'),
            new_basename => '#callgroup#.#caller#',
            add_datestamp => 1,
            suffix => ['.vcf.gz'],
          },

    };
}

sub pipeline_create_commands {
    my ($self) = @_;

    #my $sql = 'ALTER TABLE analysis_data MODIFY data MEDIUMTEXT';

    return [
        @{$self->SUPER::pipeline_create_commands},  # inheriting database and hive tables' creation
    #    $self->db_execute_command('pipeline_db', $sql),
    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},

        fai => $self->o('fai'),

        dir_label_params => ['callgroup', 'region1', 'region2'],

    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},
            '200Mb' => { 'LSF' => '-C0 -M200 -q production -R"select[mem>200] rusage[mem=200]"' },
            '500Mb' => { 'LSF' => '-C0 -M500 -q production -R"select[mem>500] rusage[mem=500]"' },
            '1Gb'   => { 'LSF' => '-C0 -M1000 -q production -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q production -R"select[mem>2000] rusage[mem=2000]"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q production -R"select[mem>4000] rusage[mem=4000]"' },
            '5Gb' => { 'LSF' => '-C0 -M5000 -q production -R"select[mem>5000] rusage[mem=5000]"' },
            '8Gb' => { 'LSF' => '-C0 -M8000 -q production -R"select[mem>8000] rusage[mem=8000]"' },
            '12Gb' => { 'LSF' => '-C0 -M12000 -q production -R"select[mem>12000] rusage[mem=12000]"' },
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
                '2->A' => { 'find_source_bams' => {'callgroup' => '#name#', 'bam_collection_id' => '#collection_id#'}},
                'A->1' => [ 'mark_seed_complete' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'find_source_bams',
            -module        => 'ReseqTrack::Hive::Process::ImportCollection',
            -meadow_type => 'LOCAL',
            -parameters    => {
                collection_id=> '#bam_collection_id#',
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
                                                'num_samples' => '#expr(scalar @{ #bam# })expr#',
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
                      1 => $self->o('call_by_samtools'),
                      2 => $self->o('call_by_gatk'),
                      3 => $self->o('call_by_freebayes'),
                      4 => 1,
                  },
                },
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -flow_into => {
                '1' => { 'call_by_samtools' => {'region2' => '#$callgroup#.#SQ_start#.#bp_start#.#SQ_end#.#bp_end#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                }},
                '2' => { 'call_by_gatk' => {'region2' => '#$callgroup#.#SQ_start#.#bp_start#.#SQ_end#.#bp_end#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                }},
                '3' => { 'call_by_freebayes' => {'region2' => '#$callgroup#.#SQ_start#.#bp_start#.#SQ_end#.#bp_end#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                }},
                '4' => [ ':////accu?bp_start=[fan_index]', ':////accu?bp_end=[fan_index]'],
            },
      });
    push(@analyses, {
            -logic_name    => 'collect_vcf',
            -module        => 'ReseqTrack::Hive::Process::BaseProcess',
            -meadow_type => 'LOCAL',
            -parameters    => {
              reseqtrack_options => {
                  flows_non_factory => {
                      1 => $self->o('call_by_samtools'),
                      2 => $self->o('call_by_gatk'),
                      3 => $self->o('call_by_freebayes'),
                      4 => 1,
                  },
                delete_param => ['bam','bai'],
              }
            },
            -flow_into => {
                1 => [ ':////accu?samtools_vcf=[fan_index]' ],
                2 => [ ':////accu?gatk_vcf=[fan_index]' ],
                3 => [ ':////accu?freebayes_vcf=[fan_index]' ],
                4 => [ ':////accu?bp_start=[fan_index]', ':////accu?bp_end=[fan_index]' ],
            },
      });
    push(@analyses, {
          -logic_name    => 'call_by_samtools',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallBySamtools',
              reference => $self->o('reference'),
              samtools => $self->o('samtools_exe'),
              bcftools => $self->o('bcftools_exe'),
              vcfutils => $self->o('vcfutils_exe'),
              bgzip => $self->o('bgzip_exe'),
              coverage_per_sample => $self->o('coverage_per_sample'),
              options => $self->o('call_by_samtools_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => { ':////accu?samtools_vcf=[fan_index]' => {'samtools_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
              -1 => [ 'call_by_samtools_himem' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_samtools_himem',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallBySamtools',
              reference => $self->o('reference'),
              samtools => $self->o('samtools_exe'),
              bcftools => $self->o('bcftools_exe'),
              vcfutils => $self->o('vcfutils_exe'),
              bgzip => $self->o('bgzip_exe'),
              coverage_per_sample => $self->o('coverage_per_sample'),
              options => $self->o('call_by_samtools_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => { ':////accu?samtools_vcf=[fan_index]' => {'samtools_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_gatk',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByGATK',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('call_by_gatk_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => { ':////accu?gatk_vcf=[fan_index]' => {'gatk_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
              -1 => [ 'call_by_gatk_himem' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_gatk_himem',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByGATK',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('call_by_gatk_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => { ':////accu?gatk_vcf=[fan_index]' => {'gatk_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_freebayes',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByFreebayes',
              reference => $self->o('reference'),
              freebayes => $self->o('freebayes_exe'),
              bgzip => $self->o('bgzip_exe'),
              options => $self->o('call_by_freebayes_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => { ':////accu?freebayes_vcf=[fan_index]' => {'freebayes_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
              -1 => [ 'call_by_freebayes_himem' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_freebayes_himem',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByFreebayes',
              reference => $self->o('reference'),
              freebayes => $self->o('freebayes_exe'),
              bgzip => $self->o('bgzip_exe'),
              options => $self->o('call_by_freebayes_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '8Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => { ':////accu?freebayes_vcf=[fan_index]' => {'freebayes_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
              -1 => [ 'call_by_freebayes_himem2' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_freebayes_himem2',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByFreebayes',
              reference => $self->o('reference'),
              freebayes => $self->o('freebayes_exe'),
              bgzip => $self->o('bgzip_exe'),
              options => $self->o('call_by_freebayes_options'),
              region_overlap => 100,
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '12Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => { ':////accu?freebayes_vcf=[fan_index]' => {'freebayes_vcf' => '#vcf#', 'fan_index' => '#fan_index#'}},
          },
      });
    push(@analyses, {
          -logic_name => 'decide_mergers',
          -module        => 'ReseqTrack::Hive::Process::BaseProcess',
          -meadow_type=> 'LOCAL',
          -parameters => {
              reseqtrack_options => {
                flows_non_factory => {
                    1 => $self->o('call_by_samtools'),
                    2 => $self->o('call_by_gatk'),
                    3 => $self->o('call_by_freebayes'),
                },
              },
          },
            -flow_into => {
                '1' => { 'merge_samtools_vcf' => {'caller' => 'samtools'}},
                '2' => { 'merge_gatk_vcf' => {'caller' => 'gatk'}},
                '3' => { 'merge_freebayes_vcf' => {'caller' => 'freebayes'}},
            },
      });


    push(@analyses, {
          -logic_name    => 'merge_samtools_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              vcf => '#samtools_vcf#',
              bgzip => $self->o('bgzip_exe'),
              reseqtrack_options => {
                decode_file_id => 'samtools_vcf',
                denestify => ['samtools_vcf','bp_start','bp_end'],
                delete_param => 'samtools_vcf',
              },
          },
          -flow_into => { '1' => [ 'store_vcf' ], },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
      });
    push(@analyses, {
          -logic_name    => 'merge_gatk_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              vcf => '#gatk_vcf#',
              bgzip => $self->o('bgzip_exe'),
              reseqtrack_options => {
                decode_file_id => 'gatk_vcf',
                denestify => ['gatk_vcf','bp_start','bp_end'],
                delete_param => 'gatk_vcf',
              },
          },
          -flow_into => { '1' => [ 'store_vcf' ], },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
      });
    push(@analyses, {
          -logic_name    => 'merge_freebayes_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              vcf => '#freebayes_vcf#',
              bgzip => $self->o('bgzip_exe'),
              reseqtrack_options => {
                decode_file_id => 'freebayes_vcf',
                denestify => ['freebayes_vcf','bp_start','bp_end'],
                delete_param => 'freebayes_vcf',
              },
          },
          -flow_into => { '1' => [ 'store_vcf' ], },
          -rc_name => '500Mb',
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

