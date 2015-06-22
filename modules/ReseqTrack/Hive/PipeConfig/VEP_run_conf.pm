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

=cut

package ReseqTrack::Hive::PipeConfig::VEP_run_conf;
use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },
        
	    'callgroup' => [],
        'pipeline_name' => 'vep',
        'final_label' => $self->o('pipeline_name'),
        'type_vcf'    => 'VCF',
        'transpose_window_size' => undef,
        'call_window_size' => undef,
        'bgzip_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
        'tabix_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/tabix',
        'vep_exe'  => '/nfs/1000g-work/G1K/work/avikd/test_hive/vep_hive/variant_effect_predictor/variant_effect_predictor.pl',
        'reference' => '/nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/bwa/grc37/human_g1k_v37.fa',
        'fai' => $self->o('reference') . '.fai',
        'vep_options' => {
                        plugin => 'AncestralAlleles,/nfs/1000g-work/G1K/work/avikd/test_hive/ancesteral_alleles_for_vep/',
                        regulatory => "",
                        offline => "",
                        vcf => "",
                        gmaf => "",
                        check_existing => "",
                        force_overwrite => '',
                        symbol => "", 
                        sift => "b",
                        polyphen => "b",
                        },
       'vep_call' => 1,
       'max_variants' => 10000,
       };
}

sub hive_meta_table {
  my ($self) = @_;
  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack' => 1,
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
            '1Gb' => { 'LSF' => '-C0 -M1000 -q production -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q production -R"select[mem>2000] rusage[mem=2000]"' },
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
                2 => [ 'find_source_vcfs' ],
            },
      });

    push(@analyses, {
            -logic_name    => 'find_source_vcfs',
            -module        => 'ReseqTrack::Hive::Process::ImportCollection',
            -meadow_type => 'LOCAL',
            -parameters    => {
                collection_type => $self->o('type_vcf'),
                collection_name=> '#callgroup#',
                output_param => 'input_vcf',
            },
            -flow_into => {
                1 => [ 'vcf_factory' ],
            },
           
      });
      push(@analyses, {
            -logic_name    => 'vcf_factory',  
            -module        => 'ReseqTrack::Hive::Process::JobFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
              factory_value => '#input_vcf#',
              temp_param_sub => { 2 => [['input_vcf', 'factory_value']]}, # temporary hack pending updates to hive code             
          },
          -flow_into => {
                2 => [ 'generate_bed' ],
            },
     });
     push(@analyses, {
          -logic_name    => 'generate_bed',
          -module        => 'ReseqTrack::Hive::Process::VcfToBed',
          -parameters    => {
              bgzip => $self->o('bgzip_exe'), 
              max_variants => $self->o('max_variants'),
              vcf =>   '#input_vcf#',
             
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
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
                max_sequences => 1,
                bed => '#vcf_bed#',
                SQ_start => '#CHROM_start#',
                SQ_end => '#CHROM_end#',

            },
		    -flow_into => {
                '2->A' => [ 'vep'],
                'A->1' => [ 'decide_mergers'],
            },
      });
#       push(@analyses, {
#             -logic_name    => 'regions_factory_2',
#             -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
#             -meadow_type => 'LOCAL',
#             -parameters    => {
#                 
#                 num_bases => $self->o('transpose_window_size'),
#                 max_sequences => 1,
#                 bed => '#vcf_bed#',
#                 SQ_start => '#SQ_start#',
#                 SQ_end => '#SQ_end#',
#             },
# 		    -flow_into => {
#                 '2' => [ 'vep',],
#                 
#             },
#       });
     push(@analyses, {
            -logic_name    => 'vep',
            -module        => 'ReseqTrack::Hive::Process::RunVep',
            -rc_name => '1Gb',
            -hive_capacity  =>  200,
            -parameters    => {
                region_overlap => 0,
                tabix_exe => $self->o('tabix_exe'),
                vcf => '#input_vcf#',
#                 bed => '#bed#',
                program_file => $self->o('vep_exe'),
                options => $self->o('vep_options'),
                temp_param_sub => { 1 => [['vep_vcf','vcf']]}, # temporary hack pending updates to hive code
            },
            -flow_into => {
              1 => [ ':////accu?vep_vcf=[fan_index]',
                      ':////accu?bp_start=[fan_index]',
                      ':////accu?bp_end=[fan_index]',
                     ],
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
                vcf => '#input_vcf#',
#                 bed => '#bed#',
                program_file => $self->o('vep_exe'),
                options => $self->o('vep_options'),
                temp_param_sub => { 1 => [['vep_vcf','vcf']]}, # temporary hack pending updates to hive code
            },
            -flow_into => {
              1 => [ ':////accu?vep_vcf=[fan_index]',
                      ':////accu?bp_start=[fan_index]',
                      ':////accu?bp_end=[fan_index]',
                     ],   
            },
            
      });
      push(@analyses, {
          -logic_name => 'decide_mergers',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              require_true => {
                  1 => $self->o('vep_call'),
                  
              },
              temp_param_sub => {
                1 => [['vcf','vep_vcf'],],
               
              }, # temporary hack pending updates to hive code
          },
            -flow_into => {
                '1' => [ 'merge_vcf' ],
                
                
            },
      });
      push(@analyses, {
          -logic_name    => 'merge_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              bgzip => $self->o('bgzip_exe'),
              
              temp_param_sub => { 1 => [['bp_start','undef'],['bp_end','undef']]}, # temporary hack pending updates to hive code
              run_tabix => 1,
              delete_param => ['vcf'],
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          
      });


return \@analyses;
}

1;
