=head1 NAME

 ReseqTrack::Hive::PipeConfig::VariantCall_conf

=head1 SYNOPSIS

  Pipeline MUST be seeded by the collection table of a ReseqTrack database (collection of bam files)
  Multiple vcf files will be created for each collection, and stored in the ReseqTrack database

  Here is an example pipeline configuration to load using reseqtrack/scripts/pipeline/load_pipeline_from_conf.pl

[variant call]
table_name=collection
seeding_module=ReseqTrack::Hive::PipeSeed::Default
config_module=ReseqTrack::Hive::PipeConfig::VariantCall_conf
config_options=-root_output_dir /path/to/dir
config_options=-reference /path/to/human.fa
  
  Options that MUST be specified in the pipeline.config_options table/column of your ReseqTrack database:

      -reference, fasta file of your reference genome.  Should be indexed for bwa and should have a .fai and .dict

  Options that have defaults but you will often want to set them in your pipeline.cofig_options table/column:

      -root_output_dir, (default is your current directory)
      -final_label, used to name your final output files (default is your pipeline name)
      -target_bed_file, for if you want to do exome calling (default undefined)
      -transpose_window_size, (default 50000000) Controls the size of the region convered by a single transposed bam
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
          depth_of_coverage => '#expr($num_samples * $coverage_per_sample)expr#',
          },

        'fai' => $self->o('reference') . '.fai',
        'target_bed_file' => undef,

        'final_label' => $self->o('pipeline_name'),
        'call_window_size' => 50000,
        'transpose_window_size' => 10000000,

        'call_by_samtools' => 1,
        'call_by_gatk' => 1,
        'call_by_freebayes' => 1,

        'vcf_type' => undef,

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
                output_columns => ['name','collection_id'],
            },
            -flow_into => {
                2 => { 'block_seed_complete' => {'callgroup' => '#name#', 'bam_collection_id' => '#collection_id#', 'ps_id' => '#ps_id#'} },
            },
      });
    push(@analyses, {
          -logic_name => 'block_seed_complete',
          -module        => 'ReseqTrack::Hive::Process::BaseProcess',
          -meadow_type=> 'LOCAL',
          -parameters => {
              flows_non_factory => [1,2],
          },
            -flow_into => {
                '2->A' => [ 'find_source_bams' ],
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
                flows_file_count_param => 'bam',
                flows_file_count => { 1 => '1+', },
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
                '2->A' => { 'transpose_bam' => {'region1' => '#expr(join(".",$callgroup,$SQ_start,$bp_start,$SQ_end,$bp_end))expr#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                'num_samples' => '#expr(scalar @{$bam})expr#',
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
                flows_factory => {
                    1 => $self->o('call_by_samtools'),
                    2 => $self->o('call_by_gatk'),
                    3 => $self->o('call_by_freebayes'),
                    4 => 1,
                },
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -flow_into => {
                '1' => { 'call_by_samtools' => {'region2' => '#expr(join(".",$callgroup,$SQ_start,$bp_start,$SQ_end,$bp_end))expr#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                }},
                '2' => { 'call_by_gatk' => {'region2' => '#expr(join(".",$callgroup,$SQ_start,$bp_start,$SQ_end,$bp_end))expr#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                }},
                '3' => { 'call_by_freebayes' => {'region2' => '#expr(join(".",$callgroup,$SQ_start,$bp_start,$SQ_end,$bp_end))expr#',
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
                flows_non_factory => {
                    1 => $self->o('call_by_samtools'),
                    2 => $self->o('call_by_gatk'),
                    3 => $self->o('call_by_freebayes'),
                    4 => 1,
                },
              delete_param => ['bam','bai'],
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
              flows_non_factory => {
                  1 => $self->o('call_by_samtools'),
                  2 => $self->o('call_by_gatk'),
                  3 => $self->o('call_by_freebayes'),
              },
          },
            -flow_into => {
                '1' => [ 'merge_samtools_vcf' ],
                '2' => [ 'merge_gatk_vcf' ],
                '3' => [ 'merge_freebayes_vcf' ],
                #'1' => { 'merge_vcf' => {'caller' => 'samtools', 'vcf' => '#samtools_vcf#'}},
                #'2' => { 'merge_vcf' => {'caller' => 'gatk', 'vcf' => '#gatk_vcf#'}},
                #'3' => { 'merge_vcf' => {'caller' => 'freebayes', 'vcf' => '#freebayes_vcf#'}},
            },
      });


    push(@analyses, {
          -logic_name    => 'merge_samtools_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              vcf => '#samtools_vcf#',
              final_label => $self->o('final_label'),
              analysis_label => '#expr(' . q('samtools.'.$final_label) . ')expr#',
              file_timestamp => 1,
              bgzip => $self->o('bgzip_exe'),
              delete_param => ['samtools_vcf'],
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
              final_label => $self->o('final_label'),
              analysis_label => '#expr(' . q('gatk.'.$final_label) . ')expr#',
              file_timestamp => 1,
              bgzip => $self->o('bgzip_exe'),
              delete_param => ['gatk_vcf'],
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
              final_label => $self->o('final_label'),
              analysis_label => '#expr(' . q('freebayes.'.$final_label) . ')expr#',
              file_timestamp => 1,
              bgzip => $self->o('bgzip_exe'),
              delete_param => ['freebayes_vcf'],
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

