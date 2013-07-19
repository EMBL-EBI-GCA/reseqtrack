=head1 NAME

 ReseqTrack::Hive::PipeConfig::VariantCall_conf

=head1 SYNOPSIS
  
  Options you MUST specify on the command line:

      -callgroup, (e.g. a population name) Can be specified multiple times.  The name of a collection of bams in your reseqtrack database
      -reference, fasta file of your reference genome.  Should be indexed for bwa and should have a .fai and .dict
      -password, for accessing the hive database
      -reseqtrack_db_name, (or -reseqtrack_db -db_name=??) your reseqtrack database

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
      -freebayes_exe, (default /nfs/1000g-work/G1K/work/bin/freebayes/bin/freebayes)

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

        callgroup => [],

        'fai' => $self->o('reference') . '.fai',
        'target_bed_file' => undef,

        'final_label' => $self->o('pipeline_name'),
        'call_window_size' => 50000,
        'transpose_window_size' => 50000000,

        'call_by_samtools' => 1,
        'call_by_gatk' => 1,
        'call_by_freebayes' => 1,


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
                  1 => $self->o('call_by_samtools'),
                  2 => $self->o('call_by_gatk'),
                  3 => $self->o('call_by_freebayes'),
              }
          },
            -flow_into => {
                '1' => [ 'call_by_samtools' ],
                '2' => [ 'call_by_gatk' ],
                '3' => [ 'call_by_freebayes' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'collect_vcf',
            -module        => 'ReseqTrack::Hive::Process::FlowDecider',
            -meadow_type => 'LOCAL',
            -parameters    => {
                require_true => {
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
              module_name => 'CallBySamtools',
              reference => $self->o('reference'),
              samtools => $self->o('samtools_exe'),
              bcftools => $self->o('bcftools_exe'),
              vcfutils => $self->o('vcfutils_exe'),
              bgzip => $self->o('bgzip_exe'),
              options => $self->o('call_by_samtools_options'),
              region_overlap => 100,
              temp_param_sub => { 1 => [['samtools_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?samtools_vcf=[fan_index]' ],
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
              options => $self->o('call_by_samtools_options'),
              region_overlap => 100,
              temp_param_sub => { 1 => [['samtools_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?samtools_vcf=[fan_index]' ],
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
              temp_param_sub => { 1 => [['gatk_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?gatk_vcf=[fan_index]' ],
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
              temp_param_sub => { 1 => [['gatk_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => [ ':////accu?gatk_vcf=[fan_index]' ],
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
              temp_param_sub => { 1 => [['freebayes_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?freebayes_vcf=[fan_index]' ],
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
              temp_param_sub => { 1 => [['freebayes_vcf','vcf']]}, # temporary hack pending updates to hive code
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  100,
          -flow_into => {
              1 => [ ':////accu?freebayes_vcf=[fan_index]' ],
          },
      });
    push(@analyses, {
          -logic_name => 'decide_mergers',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              require_true => {
                  1 => $self->o('call_by_samtools'),
                  2 => $self->o('call_by_gatk'),
                  3 => $self->o('call_by_freebayes'),
              },
              temp_param_sub => {
                1 => [['vcf','samtools_vcf'],['gatk_vcf', 'undef'],['freebayes_vcf','undef'],['caller', '"samtools"']],
                2 => [['vcf','gatk_vcf'],['samtools_vcf', 'undef'],['freebayes_vcf','undef'],['caller', '"gatk"']],
                3 => [['vcf','freebayes_vcf'],['samtools_vcf', 'undef'],['gatk_vcf','undef'],['caller', '"freebayes"']],
              }, # temporary hack pending updates to hive code
          },
            -flow_into => {
                '1' => [ 'merge_vcf' ],
                '2' => [ 'merge_vcf' ],
                '3' => [ 'merge_vcf' ],
            },
      });


    push(@analyses, {
          -logic_name    => 'merge_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              analysis_label => '#expr("call_by_".$caller)expr#',
              bgzip => $self->o('bgzip_exe'),
              delete_param => ['vcf'],
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
      });

    return \@analyses;
}

1;

