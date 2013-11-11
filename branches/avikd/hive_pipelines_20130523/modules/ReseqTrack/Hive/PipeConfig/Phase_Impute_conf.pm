=head1 NAME

 ReseqTrack::Hive::PipeConfig::Phase_Impute_conf

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
      -phase_region_bed, chrom region, required shapeit phasing, (default undefined)
      -transpose_window_size, (default 0) Controls the size of the region convered by a single transposed bam
      
      -chrom_start, chrom to start analysis (based on order in .fai file)
      -chrom_end, last chrom to include in analysis (based on order in .fai file)
      
      -chrX_string, name of chrom X in vcf/reference, i.e; “X” , “chrX” or “ChrX” (default “X”)
      -chrX_PAR1_start, starting coordinate of chrom X PAR1 region (default 1)
      -chrX_PAR1_end, ending coordinate of chrom X PAR1 region (default 2699500)
      -chrX_PAR2_start, starting coordinate of chrom X PAR2 region (default 154933000)
      -chrX_PAR2_end, ending coordinate of chrom X PAR2  region (default 155270000)
      
  Caller options
  
    BEAGLE:
      -beagle_phase , boolean, default 1, turn on/off phase by Beagle
      -gl, boolean, default 1, turn on/off considering GL from VCF file
      -phase_by_beagle_options, Beagle options 
      
    SHAPEIT:
      -shapeit_phase, boolean, default 1, turn on/off phase by Shapeit
      -phase_by_shapeit_options, Shapeit options
      -exclude_strand_flip, boolean, default 0, turn on/off excluding strand flip
      -samples_gender_info, text file with gender info for each sample, required for chrom X phasing
      
    IMPUTE:
      -impute_call, boolean, default 1, turn on/off Impute2
      -remove_haps_missing, default 1, turn on/off removing genomic windows from bed file with 0 type 2 SNPs (eg, SNPs in exome study)
      -max_base, imputation window length (default 5000000)
      -impute2_options, Impute2 options
      
      
       

  
=cut

package ReseqTrack::Hive::PipeConfig::Phase_Impute_conf;
use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },
        
	    'callgroup' => [],
        'pipeline_name' => 'phase',
        'final_label' => $self->o('pipeline_name'),
        'type_vcf'    => 'VCF',
        
        'bgzip_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
        'tabix_exe' => '/nfs/1000g-work/G1K/work/bin/tabix/tabix',
        'beagle_exe' => '/nfs/1000g-work/G1K/work/avikd/test_phase_module/b4.r1099.jar',
        'shapeit_exe' => '/nfs/1000g-work/G1K/work/avikd/PUR_run/variant_WGS/impute2_test/scripts/bin/shapeit',
        'vcfToImpute_exe'=>'/nfs/1000g-work/G1K/work/avikd/test_phase_module/vcf2impute_gen.pl',
        'impute2_exe' => "/nfs/1000g-work/G1K/work/avikd/PUR_run/variant_WGS/impute2_test/scripts/bin/impute2",
          
        'reference' => '/nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/bwa/grc37/human_g1k_v37.fa',
        'fai' => $self->o('reference') . '.fai',
        'phase_region_bed' => '/nfs/1000g-work/G1K/work/avikd/test_hive/phase_pipe/phase_run.bed',
        'samples_gender_info'=> '/nfs/1000g-work/G1K/work/avikd/test_hive/phase_pipe/export.txt',
        'reference_config' => '/nfs/1000g-work/G1K/work/avikd/test_phase_module/impute_ref_config',
        
        'chrom_start' => 22,
        'chrom_end' => 22,
        
        'transpose_window_size' => undef,
        
        'chrX_string' => 'X',
        'chrX_PAR1_start' => 1,
        'chrX_PAR1_end' => 2699500,
        'chrX_PAR2_start' => 154933000,
        'chrX_PAR2_end' => 155270000,
        

        ##beagle
        'gl' => 1,
        'beagle_phase' => 0,
        'overlap' => undef,
        'phase_by_beagle_options' => {	'window'=>24000, 
                						'overlap'=>3000, 
                						'gprobs'=>'true', },
        ##shapeit
        'exclude_strand_flip' => 0,
        'transpose_window_size' => 0,
        'shapeit_phase' => 1,
        #'phase_by_shapeit_options' => { phase => "--states 10 --window 0.5 --burn 1 --prune 1 --main 1", },
        'phase_by_shapeit_options' => { phase => "--no-mcmc", },
        
        ##impute
        'max_base' => 5000000,        
        'remove_haps_missing' => 1,
        'impute_call' => 1,
        'impute2_options' => { 'Ne'=> 20000, 'allow_large_regions' => "", },
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
         chrX_PAR1_start => $self->o('chrX_PAR1_start'),
         chrX_PAR1_end => $self->o('chrX_PAR1_end'),
         chrX_PAR2_start => $self->o('chrX_PAR2_start'),
         chrX_PAR2_start => $self->o('chrX_PAR2_start'),
         chrX_PAR2_end => $self->o('chrX_PAR2_end'),
         chrX_string => $self->o('chrX_string'),
         
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '1Gb' => { 'LSF' => '-C0 -M1000 -q production -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q production -R"select[mem>2000] rusage[mem=2000]"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q production -R"select[mem>4000] rusage[mem=4000]"' },
            '6Gb' => { 'LSF' => '-C0 -M6000 -q production -R"select[mem>6000] rusage[mem=6000]"' },
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
                output_param => 'unphased_vcf',
            },
            -flow_into => {
                1 => [ 'vcf_factory' ],
            },          
      });
      push(@analyses, {
            -logic_name    => 'vcf_factory',  
            -module        => 'ReseqTrack::Hive::Process::VcfFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
              vcf_list => '#unphased_vcf#',               
          },
          -flow_into => {
                2 => [ 'decide_callers' ],
            },
     });
    push(@analyses, {
          -logic_name => 'decide_callers',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              require_true => {
                  1 => $self->o('beagle_phase'),
                  2 => $self->o('shapeit_phase'),
              }
          },
            -flow_into => {
                '1' => [ 'regions_factory_beagle' ],
                '2' => [ 'regions_factory_shapeit' ],
            },
      });
     push(@analyses, {
            -logic_name    => 'regions_factory_beagle',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                num_bases => $self->o('transpose_window_size'),
                max_sequences => 0,
                SQ_start => $self->o('chrom_start'),
                SQ_end => $self->o('chrom_end'),                
            },
            -flow_into => {
                '2->A' => [ 'beagle_phase',],
                'A->1' => [ 'decide_vcf_mergers',],
            },       
      });

      push(@analyses, {
          -logic_name    => 'beagle_phase',
          -module        => 'ReseqTrack::Hive::Process::RunBeagle',
          -parameters    => { 
              program_file => $self->o('beagle_exe'),   
              jvm_args  => '-Xmx2g',
              beagle_options =>  $self->o('phase_by_beagle_options'), 
              gl => $self->o('gl'), 
              overlap => $self->o('overlap'),             
          },
          -rc_name => '2Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?beagle_vcf=[fan_index]',
                      ':////accu?bp_start=[fan_index]',
                      ':////accu?bp_end=[fan_index]',
                      ],
              -1 => [ 'beagle_phase_himem' ], 
         },          
      });
      push(@analyses, {
          -logic_name    => 'beagle_phase_himem',
          -module        => 'ReseqTrack::Hive::Process::RunBeagle',
          -parameters    => { 
              program_file => $self->o('beagle_exe'),   
              jvm_args  => '-Xmx4g',
              beagle_options =>  $self->o('phase_by_beagle_options'), 
              gl => $self->o('gl'),   
              overlap => $self->o('overlap'),    
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  200,
          -flow_into => {
              1 => [ ':////accu?beagle_vcf=[fan_index]',
                      ':////accu?bp_start=[fan_index]',
                      ':////accu?bp_end=[fan_index]',
                      ],
            },                
      });
      push(@analyses, {
          -logic_name => 'decide_vcf_mergers',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              require_true => {
                  1 => $self->o('beagle_phase'),
                  
              },
              temp_param_sub => {
                1 => [['vcf','beagle_vcf'],],
               
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
      push(@analyses, {
            -logic_name    => 'regions_factory_shapeit',
            -module        => 'ReseqTrack::Hive::Process::SequenceSliceFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                num_bases => $self->o('transpose_window_size'),
                max_sequences => 0,
                bed => $self->o('phase_region_bed'),
                SQ_start => $self->o('chrom_start'),
                SQ_end => $self->o('chrom_end'),                
            },
            -flow_into => {
                '2' => [ 'vcfToImpute'],

            },       
      });
      push(@analyses, {
          -logic_name    => 'vcfToImpute',
          -module        => 'ReseqTrack::Hive::Process::RunVcfToImpute',
          -parameters    => { 
              program_file => $self->o('vcfToImpute_exe'),                    
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          -flow_into => {
                '1' => [ 'modify_samples'],
            },                 
      });
      push(@analyses, {
          -logic_name    => 'modify_samples',
          -module        => 'ReseqTrack::Hive::Process::ModifyImputeSamples',
          -parameters    => { 
              samples_gender_info => $self->o('samples_gender_info'),   
              output_param => 'mod_samples',  
              delete_param => ['samples'],               
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          -flow_into => {
                '1' => [ 'shapeit_phase'],
            },                 
      });
     push(@analyses, {
          -logic_name    => 'shapeit_phase',
          -module        => 'ReseqTrack::Hive::Process::RunShapeit',
          -parameters    => { 
              program_file => $self->o('shapeit_exe'),   
              shapeit_options =>  $self->o('phase_by_shapeit_options'), 
              reference_config => $self->o('reference_config'),
              samples => '#mod_samples#',
              exclude_strand_flip => $self->o('exclude_strand_flip'),
                                     
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          -flow_into => {
                '1' => [ 'decide_impute'],
            },                 
      });
      push(@analyses, {
          -logic_name => 'decide_impute',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              require_true => {
                  1 => '#expr($impute_call && $shapeit_phase)expr#',
              },
          },
            -flow_into => {
                '1' => [ 'legends_to_bed' ],
            },
      });
      push(@analyses, {
            -logic_name    => 'legends_to_bed',
            -module        => 'ReseqTrack::Hive::Process::LegendsToBed',
            -parameters    => {
                reference_config => $self->o('reference_config'),
                max_base =>  $self->o('max_base'),
                remove_haps_missing => $self->o('remove_haps_missing'),
                delete_param => [ 'gen', 'mod_samples' ],
            },
            -rc_name => '1Gb',
            -flow_into => {
                '1' => [ 'regions_factory_2' ],
            },       
      });
      
       push(@analyses, {
            -logic_name    => 'regions_factory_2',
            -module        => 'ReseqTrack::Hive::Process::BedSliceFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                bed =>  '#legends_bed#',
            },
            -flow_into => {
                '2->A' => [ 'impute2',],
                'A->1' => [ 'decide_impute_mergers',],
            },       
      });
      push(@analyses, {
          -logic_name    => 'impute2',
          -module        => 'ReseqTrack::Hive::Process::RunImpute2',
          -parameters    => { 
              program_file => $self->o('impute2_exe'),  
              haps => '#shapeit_haps#',
              samples => '#shapeit_samples#', 
              impute2_options =>  $self->o('impute2_options'), 
              reference_config => $self->o('reference_config'),
              use_phased => 1,
              phase_result => 1,
              temp_param_sub => { 1 => [['impute2_out','impute_haps']]}, # temporary hack pending updates to hive code                           
          },
          -rc_name => '4Gb',
          -hive_capacity  =>  200,
          -flow_into => {
                1 => [ ':////accu?impute2_out=[fan_index]',
                      ':////accu?bp_start=[fan_index]',
                      ':////accu?bp_end=[fan_index]',
                     ],
                -1 => [ 'impute2_himem'],
            },
                           
      });
      push(@analyses, {
          -logic_name    => 'impute2_himem',
          -module        => 'ReseqTrack::Hive::Process::RunImpute2',
          -parameters    => { 
              program_file => $self->o('impute2_exe'),  
              haps => '#shapeit_haps#',
              samples => '#shapeit_samples#', 
              shapeit_options =>  $self->o('impute2_options'), 
              reference_config => $self->o('reference_config'),
              use_phased => 1,
              phase_result => 1,   
              temp_param_sub => { 1 => [['impute2_out','impute_haps']]}, # temporary hack pending updates to hive code                           
          },
          -rc_name => '6Gb',
          -hive_capacity  =>  200,
           -flow_into => {
                1 => [ ':////accu?impute2_out=[fan_index]',
                      ':////accu?bp_start=[fan_index]',
                      ':////accu?bp_end=[fan_index]',
                     ],
            },
                           
      });
      push(@analyses, {
          -logic_name => 'decide_impute_mergers',
          -module        => 'ReseqTrack::Hive::Process::FlowDecider',
          -meadow_type=> 'LOCAL',
          -parameters => {
              require_true => {
                  1 => $self->o('impute_call'),
                  
              },
              temp_param_sub => {
                1 => [['impute_haps','impute2_out'],],
               
              }, # temporary hack pending updates to hive code
          },
            -flow_into => {
                '1' => [ 'merge_impute' ],   
            },
      });
      push(@analyses, {
          -logic_name    => 'merge_impute',
          -module        => 'ReseqTrack::Hive::Process::MergeHaps',
          -parameters    => {
              delete_param => ['shapeit_haps', 'shapeit_samples','legends_bed', 
                               'impute', 'impute_haps', 'impute_allele_probs', 
                               'impute_info', 'impute_diplotype_ordering', 'impute_info_by_sample', 
                               'impute_samples', 'impute_summary', 'impute_warnings'],
              temp_param_sub => { 1 => [['bp_start','undef'],['bp_end','undef']]}, # temporary hack pending updates to hive code
              hap => '#impute_haps#',
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          
      });

return \@analyses;
}

1;
