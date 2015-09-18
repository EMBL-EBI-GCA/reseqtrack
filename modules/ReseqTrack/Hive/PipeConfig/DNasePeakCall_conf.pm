package ReseqTrack::Hive::PipeConfig::DNasePeakCall_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },

        'pipeline_name' => 'dnase_peak',
        seeding_module => 'ReseqTrack::Hive::PipeSeed::DNAsePeakCallSeed',
        seeding_options => {
            output_columns     => [ 'name', 'collection_id' ],
            require_columns    => $self->o( 'require_collection_columns' ),
            exclude_columns    => $self->o( 'exclude_collection_columns' ),
            require_attributes => $self->o( 'require_collection_attributes' ),
            exclude_attributes => $self->o( 'exclude_collection_attributes' ),
            metadata_file      => $self->o( 'metadata_file' ),
        },

        require_collection_columns    => {},
        exclude_collection_columns    => {},
        require_collection_attributes => {},
        exclude_collection_attributes => {},
        metadata_file                 => undef,

        hotspot_exe     => undef,
        hotspot_options => {},
        samtools        => undef,
        bedtools        => undef,
        bedToBigBedPath => undef,
        chr_file => undef,

        dir_label_params_list => [ 'sample_alias', 'experiment_source_id' ],  

        hotspot_bed_type    => undef,
        hotspot_bigbed_type => undef,
        peak_bed_type       => undef,
        collection_name     => undef,
        build_collection    => 1,


        final_output_dir     => $self->o( 'root_output_dir' ),
        final_output_layout  => undef,

        name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method => 'basic',
        hotspot_file_params => {
            new_dir       => '#final_output_dir#/#final_output_layout#',
            new_basename  => '#sample_alias#.#experiment_source_id#.Dnase.GRCh38.hotspot',
            add_datestamp => 1,
            suffix        => '.bed.gz',
        },
        peak_file_params => {
            new_dir       => '#final_output_dir#/#final_output_layout#',
            new_basename  => '#sample_alias#.#experiment_source_id#.Dnase.GRCh38.peaks',
            add_datestamp => 1,
            suffix        => '.bed.gz',
        },
       
    };
}

sub pipeline_create_commands {
    my ( $self ) = @_;

    return [
        @{ $self->SUPER::pipeline_create_commands }, 
    ];
}

sub pipeline_wide_parameters {
    my ( $self ) = @_;
    return {
        %{ $self->SUPER::pipeline_wide_parameters },


    };
}

sub resource_classes {
    my ( $self ) = @_;
    return {
            %{ $self->SUPER::resource_classes },
            '200Mb' => { 'LSF' => '-C0 -M200 -q '.$self->o( 'lsf_queue' ).' -R"select[mem>200] rusage[mem=200]"' },
            '500Mb' => { 'LSF' => '-C0 -M500 -q '.$self->o( 'lsf_queue' ).' -R"select[mem>500] rusage[mem=500]"' },
            '1Gb'   => { 'LSF' => '-C0 -M1000 -q '.$self->o( 'lsf_queue' ).' -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q '.$self->o( 'lsf_queue' ).' -R"select[mem>2000] rusage[mem=2000]"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q '.$self->o( 'lsf_queue' ).' -R"select[mem>4000] rusage[mem=4000]"' },
            '5Gb' => { 'LSF' => '-C0 -M5000 -q '.$self->o( 'lsf_queue' ).' -R"select[mem>5000] rusage[mem=5000]"' },
            '8Gb' => { 'LSF' => '-C0 -M8000 -q '.$self->o( 'lsf_queue' ).' -R"select[mem>8000] rusage[mem=8000]"' },
            '12Gb' => { 'LSF' => '-C0 -M12000 -q '.$self->o( 'lsf_queue' ).' -R"select[mem>12000] rusage[mem=12000]"' },
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
                seeding_module => $self->o( 'seeding_module' ),
                seeding_options => $self->o( 'seeding_options' ),
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
              flows_non_factory => [ 1,2 ],
            },
          },
            -flow_into => {
                '2->A' => { 'find_source_bams' => { 'callgroup' => '#name#', 'bam_collection_id' => '#collection_id#' }},
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
                1 => [ 'hotspot_peak_call' ],
            },
      });
   push(@analyses, {
            -logic_name    => 'hotspot_peak_call',
            -module        => 'ReseqTrack::Hive::Process::RunHotspot',
            -parameters    => {
                program_file => $self->o( 'hotspot_exe' ),
                options      => $self->o( 'hotspot_options' ),
            },
            -rc_name => '4Gb',
            -analysis_capacity  =>  50,
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [ 'hotspot_stats' ],
            },
      });
    push(@analyses, {
            -logic_name    => 'hotspot_stats',
            -module        => 'ReseqTrack::Hive::Process::Macs2Attribute',
            -parameters    => {
               bed      => '#hotspot_bed#',
               samtools => $self->o( 'samtools' ),
               bedtools => $self->o( 'bedtools' ),
            },
           -rc_name => '2Gb',
           -hive_capacity  =>  200,
           -flow_into => {
                1 => {  'store_bed_file' => { 'attribute_metrics' => '#attribute_metrics#',
                                              'hotspot_bed'       => '#hotspot_bed#',
                                              'peak_bed'          => '#peak_bed#',
                                            }
                     },
                },
   }); 
   push(@analyses, {
            -logic_name    => 'store_bed_file',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -parameters    => {
              type => $self->o( 'hotspot_bed_type' ),
              file => '#hotspot_bed#',
              name_file_module => $self->o( 'name_file_module' ),
              name_file_method => $self->o( 'name_file_method' ),
              name_file_params => $self->o( 'hotspot_file_params' ),
              final_output_dir => $self->o( 'final_output_dir' ),
              final_output_layout => $self->o( 'final_output_layout' ),
              collection_name => $self->o( 'collection_name' ),
              collect => $self->o( 'build_collection' ),
            },
            -rc_name => '200Mb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => {  'store_peak_file' => { 'hotspot_bed_file' => '#file#' } },
                },
   });
   push(@analyses, {
            -logic_name    => 'store_peak_file',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -parameters    => {
              type => $self->o( 'peak_bed_type' ),
              file => '#peak_bed#',
              name_file_module => $self->o( 'name_file_module' ),
              name_file_method => $self->o( 'name_file_method' ),
              name_file_params => $self->o( 'peak_file_params' ),
              final_output_dir => $self->o( 'final_output_dir' ),
              final_output_layout => $self->o( 'final_output_layout' ),
              collection_name => $self->o( 'collection_name' ),
              collect => $self->o( 'build_collection' ),
            },
            -rc_name => '200Mb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [  'bed_attributes'  ],
                },
   });
   push(@analyses, {
            -logic_name => 'bed_attributes',
            -module        => 'ReseqTrack::Hive::Process::UpdateAttribute',
            -parameters => {
                attribute_metrics =>  [ '#attribute_metrics#' ],
                collection_type => $self->o( 'hotspot_bed_type' ),
                collection_name => $self->o( 'collection_name' ),
            },
            -rc_name => '200Mb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => [  'bed_to_bigbed'  ],
                },
   });
   push(@analyses, {
            -logic_name => 'bed_to_bigbed',
            -module        => 'ReseqTrack::Hive::Process::ConvertBedToBigBed',
            -parameters => {
                bed  => [ '#hotspot_bed_file#' ],
                chr_file        => $self->o( 'chr_file' ),
                bedToBigBedPath => $self->o( 'bedToBigBedPath' ),
            },
            -rc_name => '2Gb',
            -hive_capacity  =>  200,
            -flow_into => {
                1 => {  'store_bigbed' => { 'bigbed' => '#bigbed#' } },
                },
   });
   push(@analyses, {
            -logic_name    => 'store_bigbed',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -parameters    => {
              type => $self->o( 'hotspot_bigbed_type' ),
              file => '#bigbed#',
              collection_name => $self->o( 'collection_name' ),
              collect => $self->o( 'build_collection' ),
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
