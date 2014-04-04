=head1 NAME

 ReseqTrack::Hive::PipeConfig::VEP_run_conf

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
      -type_vcf, type of vcf files to look for in the reseqtrack database, default VCF
      -final_label, used to name your final output files (default is your pipeline name)
      
     
      Caller options
      -vep_options, hash, VEP parameters
      -max_variants, max no. of variants for VEP run (default: 10000)
      
      Paths of executables:
      -bgzip_exe, (default /nfs/1000g-work/G1K/work/bin/tabix/bgzip)
      -vep_exe, (default /nfs/1000g-work/G1K/work/avikd/test_hive/vep_hive/variant_effect_predictor/variant_effect_predictor.pl)
      -tabix_exe, (default /nfs/1000g-work/G1K/work/bin/tabix/tabix)
      -vep_filter_exe, (default null)

=cut

package ReseqTrack::Hive::PipeConfig::VEP_run_conf;
use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },
        
        'pipeline_name' => 'vep',

         seeding_module => 'ReseqTrack::Hive::PipeSeed::BasePipeSeed',
         seeding_options => { 
           output_columns => ['name', 'collection_id'],
           require_columns => $self->o('require_collection_columns'),
           exclude_columns => $self->o('exclude_collection_columns'),
           require_attributes => $self->o('require_collection_attributes'),
           exclude_attributes => $self->o('exclude_collection_attributes'),
         },
  
         
        'bgzip_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
        'tabix_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/tabix',
        'vep_exe'  => '/nfs/1000g-work/G1K/work/avikd/test_dir/test_hive/vep_hive/variant_effect_predictor/variant_effect_predictor.pl',

        'vep_options' => {
                     #   plugin => 'AncestralAlleles,/nfs/1000g-work/G1K/work/avikd/test_dir/test_hive/ancesteral_alleles_for_vep/',
                        regulatory => '',
                        offline => '',
                        vcf => '',
                        gmaf => '',
                        check_existing => '',
                        force_overwrite => '',
                        symbol => '', 
                        sift => 'b',
                        polyphen => 'b',
                     #   format => 'vcf',
                        },
                        
        'vep_filter' => undef,   
        'vep_filter_exe' => undef,
        'vep_filter_options' => undef,           

        'reference' => '/nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/bwa/grc37/human_g1k_v37.fa',
        'fai' => $self->o('reference') . '.fai',
       
        'vep_call' => 1,
        'max_variants' => 10000,
        'transpose_window_size' => undef,
        
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
            new_basename => '#callgroup#.#vcf_name#.#pipeline_name#',
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
       
         dir_label_params => ['callgroup', 'vcf_name', 'region1'],
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '500Mb' => { 'LSF' => '-C0 -M500 -q production -R"select[mem>500] rusage[mem=500]"' },
            '1Gb' => { 'LSF' => '-C0 -M1000 -q production -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q production -R"select[mem>2000] rusage[mem=2000]"' },
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
               '2->A' => { 'find_source_vcfs' => {'callgroup' => '#name#', 'vcf_collection_id' => '#collection_id#'}},
               'A->1' => [ 'mark_seed_complete' ],
           },
      });
      push(@analyses, {
          -logic_name    => 'find_source_vcfs',
          -module        => 'ReseqTrack::Hive::Process::ImportCollection',
          -meadow_type => 'LOCAL',
          -parameters    => {
               collection_id=> '#vcf_collection_id#',
               output_param => 'vcf',
               reseqtrack_options => {
                  flows_file_count_param => 'vcf',
                  flows_file_count => { 1 => '1+', },
               },
          },
          -flow_into => {
                1 => [ 'vcf_factory', ':////accu?vcf=[]' ],
               
            },
      });
      push(@analyses, {
            -logic_name    => 'vcf_factory',  
            -module        => 'ReseqTrack::Hive::Process::JobFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
              factory_value => '#vcf#',                
          },
          -flow_into => {
                2 => { 'generate_bed' => { 'vcf' => '#factory_value#', 'vcf_name' => '#expr( ($factory_value =~ /([^\/]*)\.vcf(\.gz)?$/)[0] )expr#' }},
            },
     });
     push(@analyses, {
          -logic_name    => 'generate_bed',
          -module        => 'ReseqTrack::Hive::Process::VcfToBed',
          -parameters    => {
              bgzip => $self->o('bgzip_exe'), 
              max_variants => $self->o('max_variants'),             
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          -flow_into => {  
                '1' => { 'regions_factory_1' => { 'vcf_bed' => '#bed#', 'vcf_name' => '#vcf_name#' ,'CHROM_start' => '#CHROM_start#', 'CHROM_end' => '#CHROM_end#' , },},
            },
      });
      push(@analyses, {
            -logic_name    => 'regions_factory_1',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                num_bases => $self->o('transpose_window_size'),
                max_sequences => 1,
                bed => '#vcf_bed#',
                SQ_start => '#CHROM_start#',
                SQ_end => '#CHROM_end#',
              },  
            -flow_into => {
                 '2->A' =>  { 'vep' => { 'region1' => '#callgroup#.#vcf_name#.#SQ_start#.#bp_start#.#SQ_end#.#bp_end#',
                              'SQ_start' => '#SQ_start#', 'bp_start' => '#bp_start#', 
                              'SQ_end' => '#SQ_end#', 'bp_end' => '#bp_end#','fan_index' => '#fan_index#' },
                             },
                 'A->1' =>  [ 'decide_mergers' ],
            },
      });
     push(@analyses, {
            -logic_name    => 'vep',
            -module        => 'ReseqTrack::Hive::Process::RunVep',
            -rc_name => '1Gb',
            -hive_capacity  =>  200,
            -parameters    => {
                region_overlap => 0,
                tabix_exe => $self->o('tabix_exe'),
                program_file => $self->o('vep_exe'),
                options => $self->o('vep_options'),
                vep_filter => $self->o('vep_filter'),
                vep_filter_exe => $self->o('vep_filter_exe'),
                vep_filter_options => $self->o('vep_filter_options'),
                
                reseqtrack_options => {
                  encode_file_id => 'vcf',
                },
            },
            -flow_into => {
              1 => { ':////accu?vcf=[fan_index]' => {'vcf' => '#vcf#', 'fan_index' => '#fan_index#'},
                      ':////accu?bp_start=[fan_index]' => {'bp_start' => '#bp_start#', 'fan_index' => '#fan_index#'},
                      ':////accu?bp_end=[fan_index]' => {'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#'},},                
              -1 => [ 'vep_himem' ],
            },
            
      });
      push(@analyses, {
            -logic_name    => 'vep_himem',
            -module        => 'ReseqTrack::Hive::Process::RunVep',
            -rc_name => '2Gb',
            -hive_capacity  =>  200,
            -parameters    => {
                region_overlap => 0,
                tabix_exe => $self->o('tabix_exe'),
                program_file => $self->o('vep_exe'),
                options => $self->o('vep_options'),
                vep_filter => $self->o('vep_filter'),
                vep_filter_exe => $self->o('vep_filter_exe'),
                vep_filter_options => $self->o('vep_filter_options'),
                reseqtrack_options => {
                  encode_file_id => 'vcf',
                },
            },
            -flow_into => {
              1 => { ':////accu?vcf=[fan_index]' => {'vcf' => '#vcf#', 'fan_index' => '#fan_index#'},
                      ':////accu?bp_start=[fan_index]' => {'bp_start' => '#bp_start#', 'fan_index' => '#fan_index#'},
                      ':////accu?bp_end=[fan_index]' => {'bp_end' => '#bp_end#', 'fan_index' => '#fan_index#'},},
            },   
      });
      push(@analyses, {
          -logic_name => 'decide_mergers',
          -module        => 'ReseqTrack::Hive::Process::BaseProcess',
          -meadow_type=> 'LOCAL',
          -parameters => {
              reseqtrack_options => {
                flows_non_factory => {
                    1 => $self->o('vep_call'),
                },
              },
          },
            -flow_into => {
                '1' => { 'merge_vcf'=> {'caller' => 'vep'} },  
            },
      });
      push(@analyses, {
          -logic_name    => 'merge_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              vcf => '#vcf#',
              bgzip => $self->o('bgzip_exe'),
              reseqtrack_options => {
                decode_file_id => 'vcf',
                #denestify => ['vep_vcf','bp_start','bp_end'],
                delete_param => ['vcf','vcf_bed'],
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
          -rc_name => '500Mb',
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
