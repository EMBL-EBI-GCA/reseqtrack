package AccessibleGenome::BamUtil_conf;

use strict;
use ReseqTrack::Tools::Exception;
use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },

        'pipeline_name' => 'accessible_genome',

        seeding_module => 'AccessibleGenome::BamUtilPipeSeed',
        seeding_options => {
            output_columns => $self->o('output_columns'),
            require_columns => $self->o('require_file_table_columns')
		},
    	
        'output_columns' =>['name'],
        'require_file_table_columns' => {},

        'bgzip_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
        'bam_util_exe' => '/nfs/1000g-work/G1K/work/bin/bamUtil/bin/bam',

        bam_util_options => {}, 
        
        final_output_dir => $self->o('final_output_dir'),
        final_output_layout => '#POP#',
        name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method => 'basic',
        name_file_params => {
            new_dir => '#final_output_dir#/#final_output_layout#',
            new_basename => '#file_basename#.stats',
            add_datestamp => 0,
            suffix => '.gz',
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
        %{$self->SUPER::pipeline_wide_parameters},

        dir_label_params => ['POP', 'file_basename', 'region1', 'region2'],

    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},

            '200Mb' => { 'LSF' => '-C0 -M200 -q '.$self->o('lsf_queue').' -R"select[mem>200] rusage[mem=200] select[' .$self->o('lsf_resource') . ']"' },
            '500Mb' => { 'LSF' => '-C0 -M500 -q '.$self->o('lsf_queue').' -R"select[mem>500] rusage[mem=500] select[' .$self->o('lsf_resource') . ']"' },
            '800Mb' => { 'LSF' => '-C0 -M800 -q '.$self->o('lsf_queue').' -R"select[mem>800] rusage[mem=800] select[' .$self->o('lsf_resource') . ']"' },
            '1Gb'   => { 'LSF' => '-C0 -M1000 -q '.$self->o('lsf_queue').' -R"select[mem>1000] rusage[mem=1000] select[' .$self->o('lsf_resource') . ']"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q '.$self->o('lsf_queue').' -R"select[mem>2000] rusage[mem=2000] select[' .$self->o('lsf_resource') . ']"' },
            '3Gb' => { 'LSF' => '-C0 -M3000 -q '.$self->o('lsf_queue').' -R"select[mem>3000] rusage[mem=3000] select[' .$self->o('lsf_resource') . ']"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q '.$self->o('lsf_queue').' -R"select[mem>4000] rusage[mem=4000] select[' .$self->o('lsf_resource') . ']"' },
            '5Gb' => { 'LSF' => '-C0 -M5000 -q '.$self->o('lsf_queue').' -R"select[mem>5000] rusage[mem=5000] select[' .$self->o('lsf_resource') . ']"' },
            '6Gb' => { 'LSF' => '-C0 -M6000 -q '.$self->o('lsf_queue').' -R"select[mem>6000] rusage[mem=6000] select[' .$self->o('lsf_resource') . ']"' },
            '8Gb' => { 'LSF' => '-C0 -M8000 -q '.$self->o('lsf_queue').' -R"select[mem>8000] rusage[mem=8000] select[' .$self->o('lsf_resource') . ']"' },
            '12Gb' => { 'LSF' => '-C0 -M12000 -q '.$self->o('lsf_queue').' -R"select[mem>12000] rusage[mem=12000] select[' .$self->o('lsf_resource') . ']"' },
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
                2 => [ 'regions_factory_1' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'regions_factory_1',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -parameters    => {
                num_bases => 100000000,
                max_sequences => 200,
            },
            -flow_into => {
                '2->A' => { 'regions_factory_2' => {'region1' => '#file_basename#.#SQ_start#.#bp_start#.#SQ_end#.#bp_end#',
                                                'SQ_start' => '#SQ_start#','SQ_end' => '#SQ_end#','bp_start' => '#bp_start#','bp_end'=>'#bp_end#','fan_index' => '#fan_index#',
                                                }},
                'A->1' => [ 'merge_stats'],
            },
      });
    push(@analyses, {
            -logic_name    => 'regions_factory_2',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -parameters    => {
                num_bases => 5000000,
                max_sequences => 1,
            },
            -rc_name => '200Mb',
            -analysis_capacity  =>  4,
            -hive_capacity  =>  200,
            -flow_into => {
                '2->A' => { 'bam_util' => {'region2' => '#file_basename#.#SQ_start#.#bp_start#.#bp_end#',
                                                'SQ' => '#SQ_start#', 'bp_start' => '#bp_start#', 'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#',
                                                },
                },
                'A->1' => [ 'collect_stats'],
            },
      });
    push(@analyses, {
            -logic_name    => 'collect_stats',
            -module        => 'ReseqTrack::Hive::Process::BaseProcess',
            -meadow_type => 'LOCAL',
            -flow_into => {
                1 => [ ':////accu?stats=[fan_index]' ],
            },
      });
    push(@analyses, {
          -logic_name    => 'bam_util',
          -module        => 'AccessibleGenome::BamUtilHiveProcess',
          -parameters    => {
              bam_util_exe => $self->o('bam_util_exe'),
              options => $self->o('bam_util_options'),
              reseqtrack_options => {
                encode_file_id => 'stats',
              },
          },
          -rc_name => '6Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?stats=[fan_index]' ],
              -1 => [ 'bam_util_himem'],
          },
      });
    push(@analyses, {
          -logic_name    => 'bam_util_himem',
          -module        => 'AccessibleGenome::BamUtilHiveProcess',
          -parameters    => {
              bam_util_exe => $self->o('bam_util_exe'),
              options => $self->o('bam_util_options'),
              reseqtrack_options => {
                encode_file_id => 'stats',
              },
          },
          -rc_name => '8Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?stats=[fan_index]' ],
          },
      });


    push(@analyses, {
          -logic_name    => 'merge_stats',
          -module        => 'AccessibleGenome::MergeStats',
          -parameters    => {
              bgzip => $self->o('bgzip_exe'),
              run_tabix => 1,
              reseqtrack_options => {
                decode_file_id => 'stats',
                denestify => 'stats',
                delete_param => 'stats',
              }
          },
          -flow_into => { '1' => [ 'move_file' ], },
          -rc_name => '500Mb',
          -hive_capacity  =>  200,
      });
    push(@analyses, {
            -logic_name    => 'move_file',
            -module        => 'AccessibleGenome::MoveFile',
            -parameters    => {
              file => '#stats#',
              name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
              name_file_method => 'basic',
              name_file_params => $self->o('name_file_params'),
              final_output_dir => $self->o('final_output_dir'),
              final_output_layout => $self->o('final_output_layout'),
            },
            -meadow_type => 'LOCAL',
            -flow_into => {1 => {'move_tbi' => {'stats' => '#file#'}}},
      });
    push(@analyses, {
            -logic_name    => 'move_tbi',
            -module        => 'AccessibleGenome::MoveFile',
            -parameters    => {
              file => '#tbi#',
              name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
              name_file_method => 'basic',
              name_file_params => {new_full_path => '#stats#.tbi'},
            },
            -meadow_type => 'LOCAL',
            -flow_into => {1 => {'mark_seed_complete' => {'tbi' => '#file#'}}},
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

