=head1 NAME

 ReseqTrack::Hive::PipeConfig::LobSTRAllelotype_conf

=head1 SYNOPSIS

  Pipeline must be seeded by the collection table of a ReseqTrack database (collection of lobstr bam files)
  A vcf file will be created for each collection, and stored in the ReseqTrack database

  Here is an example pipeline configuration to load using reseqtrack/scripts/pipeline/load_pipeline_from_conf.pl

[lobstr_allelotype]
table_name=collection
config_module=ReseqTrack::Hive::PipeConfig::LobSTRAllelotype_conf
config_options=-root_output_dir /path/to/dir
config_options=-reference /path/to/human.fa
config_options=-callgroup_type LOBSTR_BAMS
config_options=-lobstr_ref_index_prefix /path/to/resources/lobstr_
config_options=-lobstr_str_info /path/to/resources/lobstr_v2.0.3_hg19_strinfo.tab
config_options=-lobstr_noise_model /path/to/models/illumina_v2.0.3
  
  Options that MUST be specified in the pipeline.config_options table/column of your ReseqTrack database:

      -callgroup_type, type of collection (of bams) in the reseqtrack database used for seeding the pipeline
      -reference, fasta file of your reference genome.  Should have a .fai
      -lobstr_ref_index_prefix, e.g. /path/to/resources/lobstr_
      -lobstr_str_info, e.g. /path/to/resources/lobstr_v2.0.3_hg19_strinfo.tab
      -lobstr_noise_model, e.g. /path/to//models/illumina_v2.0.3

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

      -call_by_lobstr_options, key/val pairs, passed on to the CallByLobstr module

      Paths of executables:
      -lobstr_exe, (default /nfs/1000g-work/G1K/work/bin/lobstr/bin/allelotype)
      -bgzip_exe, (default /nfs/1000g-work/G1K/work/bin/tabix/bgzip)

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
            '200Mb' => { 'LSF' => '-C0 -M200 -q '.$self->o('lsf_queue').' -R"select[mem>200] rusage[mem=200]"' },
            '500Mb' => { 'LSF' => '-C0 -M500 -q '.$self->o('lsf_queue').' -R"select[mem>500] rusage[mem=500]"' },
            '1Gb'   => { 'LSF' => '-C0 -M1000 -q '.$self->o('lsf_queue').' -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q '.$self->o('lsf_queue').' -R"select[mem>2000] rusage[mem=2000]"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q '.$self->o('lsf_queue').' -R"select[mem>4000] rusage[mem=4000]"' },
            '5Gb' => { 'LSF' => '-C0 -M5000 -q '.$self->o('lsf_queue').' -R"select[mem>5000] rusage[mem=5000]"' },
            '8Gb' => { 'LSF' => '-C0 -M8000 -q '.$self->o('lsf_queue').' -R"select[mem>8000] rusage[mem=8000]"' },
            '12Gb' => { 'LSF' => '-C0 -M12000 -q '.$self->o('lsf_queue').' -R"select[mem>12000] rusage[mem=12000]"' },
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
                '2->A' => { 'call_by_lobstr' => {'region2' => '#callgroup#.#SQ_start#',
                                                'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                },
                },
                'A->1' => [ 'collect_vcf'],
            },
      });
    push(@analyses, {
            -logic_name    => 'collect_vcf',
            -module        => 'ReseqTrack::Hive::Process::BaseProcess',
            -meadow_type => 'LOCAL',
            -parameters    => {
              reseqtrack_options => {
                delete_param => ['bam','bai'],
              },
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
              reseqtrack_options => {
                encode_file_id => 'vcf',
              },
          },
          -rc_name => '500Mb',
          -hive_capacity  =>  50,
          -flow_into => {
              1 => [ ':////accu?vcf=[fan_index]', ':////accu?bp_start=[fan_index]', ':////accu?bp_end=[fan_index]', ],
          },
      });

    push(@analyses, {
          -logic_name    => 'merge_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              bgzip => $self->o('bgzip_exe'),
              reseqtrack_options => {
                decode_file_id => 'vcf',
                denestify => ['vcf','bp_start','bp_end'],
                delete_param => 'vcf',
              }
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

