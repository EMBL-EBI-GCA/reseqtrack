=head1 NAME

 ReseqTrack::Hive::PipeConfig::VQSR_conf

=head1 SYNOPSIS
  
  Options you MUST specify on the command line:

      -callgroup, (e.g. a population name) Can be specified multiple times.  The name of a collection of bams in your reseqtrack database
      -reference, fasta file of your reference genome.  Should be indexed for bwa and should have a .fai and .dict
      -password, for accessing the hive database
      -reseqtrack_db_name, (or -reseqtrack_db -db_name=??) your reseqtrack database

      -snp_alleles, path to a vcf file containing sites for recalibration (should be tabixed)
      -indel_alleles, path to a vcf file containing sites for recalibration (should be tabixed)

      -snp_resources, key/val pairs e.g. omni='known=true,training=true,truth=true,prior=12.0 /path/to/file.vcf
      -indel_resources, key/val pairs e.g.  mills='known=true,training=true,truth=true,prior=12.0 /path/to/file.vcf
          Resource files should be tabixed

  Options that have defaults but you will often want to modify:

      Connection to the hive database:
      -pipeline_db -host=???, (default mysql-g1k)
      -pipeline_db -port=???, (default 4175)
      -pipeline_db -user=???, must have write access (default g1krw)
      -dipeline_db -dbname=???, (default is a mixture of your unix user name + the pipeline name)

      Connection to the reseqtrack database:
      -reseqtrack_db -host=???, (default mysql-g1k)
      -reseqtrack_db -user=???, read only access is OK (default g1kro)
      -reseqtrack_db -port=???, (default 4175)
      -reseqtrack_db -pass=???, (default undefined)

      -root_output_dir, (default is your current directory)
      -type_bam, type of bam files to look for in the reseqtrack database, default BAM
      -final_label, used to name your final output files (default is your pipeline name)
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

        'type_bam'    => 'BAM',
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

        callgroup => [],

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

        'final_label' => $self->o('pipeline_name'),
        'call_window_size' => 50000,
        'transpose_window_size' => 50000000,

        'vqsr_snps' => 1,
        'vqsr_indels' => 1,


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
    };
}

sub pipeline_analyses {
    my ($self) = @_;

    my @analyses;
    push(@analyses, {
            -logic_name    => 'callgroups_factory',
            -module        => 'ReseqTrack::Hive::Process::JobFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                factory_value => '#callgroup#',
                temp_param_sub => { 2 => [['callgroup','factory_value']]}, # temporary hack pending updates to hive code
            },
            -input_ids => [{callgroup => $self->o('callgroup')}],
            -flow_into => {
                2 => [ 'find_source_bams' ],
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
                temp_param_sub => { 1 => [['bam','undef']]}, # temporary hack pending updates to hive code
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
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -flow_into => {
                '2' => [ 'decide_callers',
                          ':////accu?bp_start=[fan_index]',
                          ':////accu?bp_end=[fan_index]',
                          ],
            },
      });
    push(@analyses, {
          -logic_name => 'decide_callers',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              require_true => {
                  1 => $self->o('vqsr_snps'),
                  2 => $self->o('vqsr_indels'),
              }
          },
            -flow_into => {
                '1' => [ 'call_by_gatk_snps' ],
                '2' => [ 'call_by_gatk_indels' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'collect_vcf',
            -module        => 'ReseqTrack::Hive::Process::FlowDecider',
            -meadow_type => 'LOCAL',
            -parameters    => {
                require_true => {
                  1 => $self->o('vqsr_snps'),
                  2 => $self->o('vqsr_indels'),
                  3 => 1,
              },
              delete_param => ['bam','bai'],
            },
            -flow_into => {
                1 => [ ':////accu?gatk_snps_vcf=[fan_index]'],
                2 => [ ':////accu?gatk_indels_vcf=[fan_index]'],
                3 => [ ':////accu?bp_start=[fan_index]',
                       ':////accu?bp_end=[fan_index]',
                    ],
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
              temp_param_sub => { 1 => [['gatk_snps_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?gatk_snps_vcf=[fan_index]' ],
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
              temp_param_sub => { 1 => [['gatk_snps_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => [ ':////accu?gatk_snps_vcf=[fan_index]' ],
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
              temp_param_sub => { 1 => [['gatk_indels_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?gatk_indels_vcf=[fan_index]' ],
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
              extra_options => $self->o('call_by_gatk_indels_options'),
              alleles_file => $self->o('indel_alleles'),
              options => '#expr({%$extra_options, glm => "INDEL", gt_mode => "GENOTYPE_GIVEN_ALLELES", out_mode => "EMIT_ALL_SITES", alleles => $alleles_file})expr#',
              region_overlap => 100,
              temp_param_sub => { 1 => [['gatk_indels_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => [ ':////accu?gatk_indels_vcf=[fan_index]' ],
          },
      });
    push(@analyses, {
          -logic_name => 'decide_mergers',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              require_true => {
                  1 => $self->o('vqsr_snps'),
                  2 => $self->o('vqsr_indels'),
              },
              temp_param_sub => {
                1 => [['vcf','gatk_snps_vcf'],['gatk_indels_vcf', 'undef']],
                2 => [['vcf','gatk_indels_vcf'],['gatk_snps_vcf', 'undef']],
              }, # temporary hack pending updates to hive code
          },
            -flow_into => {
                '1' => [ 'merge_vcf_snps' ],
                '2' => [ 'merge_vcf_indels' ],
            },
      });


    push(@analyses, {
          -logic_name    => 'merge_vcf_snps',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              bgzip => $self->o('bgzip_exe'),
              delete_param => ['vcf'],
              temp_param_sub => { 1 => [['bp_start','undef'],['bp_end','undef']]}, # temporary hack pending updates to hive code
              run_tabix => 1,
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
              delete_param => ['vcf'],
              temp_param_sub => { 1 => [['bp_start','undef'],['bp_end','undef']]}, # temporary hack pending updates to hive code
              run_tabix => 1,
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
              mode => 'SNP',
          },
          -rc_name => '2Gb',
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
              mode => 'INDEL',
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ 'apply_recalibration_indels' ],
          },
      });

    push(@analyses, {
          -logic_name    => 'apply_recalibration_snps',
          -module        => 'ReseqTrack::Hive::Process::RunApplyRecalibration',
          -parameters    => {
              delete_param => ['vcf', 'recal_file', 'tbi'],
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('apply_recalibration_snps_options'),
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
      });
    push(@analyses, {
          -logic_name    => 'apply_recalibration_indels',
          -module        => 'ReseqTrack::Hive::Process::RunApplyRecalibration',
          -parameters    => {
              delete_param => ['vcf', 'recal_file', 'tbi'],
              reference => $self->o('reference'),
              gatk_dir => $self->o('gatk_dir'),
              options => $self->o('apply_recalibration_indels_options'),
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
      });

    return \@analyses;
}

1;

