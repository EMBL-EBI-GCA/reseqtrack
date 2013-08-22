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
      -type_bam, type of bam files to look for in the reseqtrack database, default BAM
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
        'target_bed_file' => undef,

        'final_label' => $self->o('pipeline_name'),
        'call_window_size' => 50000,
        'transpose_window_size' => 50000000,

        'call_by_samtools' => 1,
        'call_by_gatk' => 1,
        'call_by_freebayes' => 1,

        'vcf_type' => undef,

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
        %{$self->SUPER::pipeline_wide_parameters},

        fai => $self->o('fai'),

        labels => ['#callgroup#',
                   '#expr($region1 ? $callgroup.\'.\'.$region1 : undef)expr',
                   '#expr($region2 ? $callgroup.\'.\'.$region2 : undef)expr',
                   ],

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
                output_columns => 'name',
            },
            -flow_into => {
                2 => { 'block_seed_complete' => {'callgroup' => '#name#'} },
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
                collection_type => $self->o('type_bam'),
                collection_name=> '#callgroup#',
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
                '2->A' => [ 'transpose_bam' ],
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
                region1 => '#expr(join(".",$callgroup,$SQ_start.$bp_start,$SQ_end,$bp_end))expr#)'
            },
            -rc_name => '2Gb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -flow_into => {
                '1->A' => { 'regions_factory_2' => {'region1' => '#region1#'}},
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
                '1' => [ 'call_by_samtools' ],
                '2' => [ 'call_by_gatk' ],
                '3' => [ 'call_by_freebayes' ],
                '4' => [ ':////accu?bp_start=[fan_index]', ':////accu?bp_end=[fan_index]', ],
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
                1 => [ ':////accu?samtools_vcf=[fan_index]'],
                2 => [ ':////accu?gatk_vcf=[fan_index]'],
                3 => [ ':////accu?freebayes_vcf=[fan_index]'],
                4 => [ ':////accu?bp_start=[fan_index]',
                       ':////accu?bp_end=[fan_index]',
                    ],
            },
      });
    push(@analyses, {
          -logic_name    => 'call_by_samtools',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              region2 => '#expr(join(".",$callgroup,$SQ_start.$bp_start,$SQ_end,$bp_end))expr#)',
              module_name => 'CallBySamtools',
              reference => $self->o('reference'),
              samtools => $self->o('samtools_exe'),
              bcftools => $self->o('bcftools_exe'),
              vcfutils => $self->o('vcfutils_exe'),
              bgzip => $self->o('bgzip_exe'),
              options => $self->o('call_by_samtools_options'),
              region_overlap => 100,
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => { ':////accu?samtools_vcf=[fan_index]' => {'samtools_vcf' => '#vcf#'}},
              -1 => [ 'call_by_samtools_himem' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_samtools_himem',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              region2 => '#expr(join(".",$callgroup,$SQ_start.$bp_start,$SQ_end,$bp_end))expr#)',
              module_name => 'CallBySamtools',
              reference => $self->o('reference'),
              samtools => $self->o('samtools_exe'),
              bcftools => $self->o('bcftools_exe'),
              vcfutils => $self->o('vcfutils_exe'),
              bgzip => $self->o('bgzip_exe'),
              options => $self->o('call_by_samtools_options'),
              region_overlap => 100,
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => { ':////accu?samtools_vcf=[fan_index]' => {'samtools_vcf' => '#vcf#'}},
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_gatk',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              region2 => '#expr(join(".",$callgroup,$SQ_start.$bp_start,$SQ_end,$bp_end))expr#)',
              module_name => 'CallByGATK',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('call_by_gatk_options'),
              region_overlap => 100,
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => { ':////accu?gatk_vcf=[fan_index]' => {'gatk_vcf' => '#vcf#'}},
              -1 => [ 'call_by_gatk_himem' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_gatk_himem',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              region2 => '#expr(join(".",$callgroup,$SQ_start.$bp_start,$SQ_end,$bp_end))expr#)',
              module_name => 'CallByGATK',
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('call_by_gatk_options'),
              region_overlap => 100,
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => { ':////accu?gatk_vcf=[fan_index]' => {'gatk_vcf' => '#vcf#'}},
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_freebayes',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              region2 => '#expr(join(".",$callgroup,$SQ_start.$bp_start,$SQ_end,$bp_end))expr#)',
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
              1 => { ':////accu?freebayes_vcf=[fan_index]' => {'freebayes_vcf' => '#vcf#'}},
              -1 => [ 'call_by_freebayes_himem' ],
          },
      });
    push(@analyses, {
          -logic_name    => 'call_by_freebayes_himem',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              region2 => '#expr(join(".",$callgroup,$SQ_start.$bp_start,$SQ_end,$bp_end))expr#)',
              module_name => 'CallByFreebayes',
              reference => $self->o('reference'),
              freebayes => $self->o('freebayes_exe'),
              bgzip => $self->o('bgzip_exe'),
              options => $self->o('call_by_freebayes_options'),
              region_overlap => 100,
              temp_param_sub => { 1 => [['freebayes_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => { ':////accu?freebayes_vcf=[fan_index]' => {'freebayes_vcf' => '#vcf#'}},
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
                '1' => { 'merge_vcf' => {'caller' => 'samtools', 'vcf' => '#samtools_vcf#'}},
                '2' => { 'merge_vcf' => {'caller' => 'gatk', 'vcf' => '#gatk_vcf#'}},
                '3' => { 'merge_vcf' => {'caller' => 'freebayes', 'vcf' => '#freebayes_vcf#'}},
            },
      });


    push(@analyses, {
          -logic_name    => 'merge_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              final_label => $self->o('final_label'),
              analysis_label => '#expr(' . q($caller.'.'.$final_label) . ')expr#',
              bgzip => $self->o('bgzip_exe'),
              delete_param => ['vcf'],
          },
          -flow_into => { '1' => [ 'store_vcf' ], },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
      });
    push(@analyses, {
            -logic_name    => 'store_vcf',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -meadow_type => 'LOCAL',
            -parameters    => {
              type => $self->o('vcf_type'),
              file => '#vcf#',
            },
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

