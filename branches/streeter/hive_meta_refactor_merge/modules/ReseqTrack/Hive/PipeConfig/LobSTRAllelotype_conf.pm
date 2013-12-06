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



package ReseqTrack::Hive::PipeConfig::LobSTRAllelotype_conf;

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

        'bgzip_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
        'lobstr_exe' => '/nfs/1000g-work/G1K/work/bin/lobstr/bin/allelotype',

        call_by_lobstr_options => {}, # use module defaults

        'fai' => $self->o('reference') . '.fai',
        'target_bed_file' => undef,

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
            new_basename => '#callgroup#.lobstr',
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
              flows_non_factory => [1,2],
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
                num_bases => 99999999999,
                max_sequences => 200,
                bed => $self->o('target_bed_file'),
            },
            -flow_into => {
                '2->A' => { 'regions_factory_2' => {'region1' => '#callgroup#.#SQ_start#.#SQ_end#',
                                                'SQ_start' => '#SQ_start#','SQ_end' => '#SQ_end#','fan_index' => '#fan_index#',
                                                }},
                'A->1' => [ 'merge_vcf'],
            },
      });
    push(@analyses, {
            -logic_name    => 'regions_factory_2',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -parameters    => {
                num_bases => 0,
                max_sequences => 1,
                bed => $self->o('target_bed_file'),
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -flow_into => {
                '2' => { 'call_by_lobstr' => {'region2' => '#callgroup#.#SQ_start#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                },
                          ':////accu?bp_start=[fan_index]' => {bp_start => '#bp_start#', 'fan_index' => '#fan_index#'},
                          ':////accu?bp_end=[fan_index]' => {bp_end => '#bp_end#', 'fan_index' => '#fan_index#'},
                }
            },
      });
    push(@analyses, {
            -logic_name    => 'collect_vcf',
            -module        => 'ReseqTrack::Hive::Process::BaseProcess',
            -meadow_type => 'LOCAL',
            -parameters    => {
              delete_param => ['bam','bai'],
            },
            -flow_into => {
                1 => [ ':////accu?vcf=[fan_index]', ':////accu?bp_start=[fan_index]', ':////accu?bp_end=[fan_index]' ],
            },
      });
    push(@analyses, {
          -logic_name    => 'call_by_lobstr',
          -module        => 'ReseqTrack::Hive::Process::RunVariantCall',
          -parameters    => {
              module_name => 'CallByLobSTR',
              ref_index_prefix => $self->o('lobstr_ref_index_prefix'),
              lobstr => $self->o('lobstr_exe'),
              noise_model => $self->o('lobstr_noise_model'),
              str_info => $self->o('lobstr_str_info'),
              options => $self->o('call_by_lobstr_options'),
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?vcf=[fan_index]' ],
          },
      });

    push(@analyses, {
          -logic_name    => 'merge_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
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

