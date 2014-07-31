=head1 NAME

 ReseqTrack::Hive::PipeConfig::AnnoVcf_conf

=head1 SYNOPSIS

  Pipeline must be seeded by the collection table of a ReseqTrack database (collection of genotype VCF files, one chromosome per collection)
  A AnnoVcf matrix .txt file will be created for each collection (chromosome), containing coverage info for all samples, and stored in the ReseqTrack database

  Here is an example pipeline configuration to load using reseqtrack/scripts/pipeline/load_pipeline_from_conf.pl

[AnnoVcf]
table_name=file
config_options=-file_type VCF
config_options=-bam_type EXOME
config_module=ReseqTrack::Hive::PipeConfig::AnnoVcf_conf
config_options=-root_output_dir /path/to/dir
  
  Options that MUST be specified in the pipeline.config_options table/column of your ReseqTrack database:

      -file_type, type of file in the reseqtrack database used for seeding the pipeline

  Options that have defaults but you will often want to set them in your pipeline.cofig_options table/column:

      -seeding_module, (default is ReseqTrack::Hive::PipeSeed::BasePipeSeed) override this with a project-specific module
      -seeding_options, hashref passed to the seeding module.  Override the defaults only if using a different seeding module.

      -require_file_columns, -exclude_file_columns, -require_file_attributes, -exclude_file_attributes'
            These are hashrefs, add to these to control which files are used to seed the pipeline
            e.g. -exclude_file_columns name=GBR

      -name_file_module, (default is ReseqTrack::Hive::NameFile::BaseNameFile) override this with a project-specific module. Controls how your output bam file is named.
      -name_file_method, (default is 'basic'), controls which subroutine of your name_file_module is used to name your bam file.
      -final_output_dir, (default is your root_output_dir) the root output directory for your final bam file
      -name_file_params. a hash ref, passed to your name_file_module.  Change the default values to control how your final output file is named.

      -root_output_dir, (default is your current directory) This is where working files go, i.e. not necessarily the final resting place of your output vcf
      

      Paths of executables:
      -tabix, (default is /nfs/1000g-work/G1K/work/bin/tabix/tabix)
      -BedTools_program, (default is /nfs/1000g-work/G1K/work/bin/BEDTools/bin/multiBamCov)
	  -paste_exe (default is /usr/bin/paste)

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


package ReseqTrack::Hive::PipeConfig::AnnoVcf_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },

        'pipeline_name' => 'AnnoVcf',

        seeding_module => 'ReseqTrack::Hive::PipeSeed::BasePipeSeed',
        seeding_options => {
            output_columns => ['name', 'file_id'], 
            require_columns => $self->o('require_file_columns'),
            exclude_columns => $self->o('exclude_file_columns'),
            require_attributes => $self->o('require_file_attributes'),
            exclude_attributes => $self->o('exclude_file_attributes'),
          },
        		
        'require_file_columns' => {'type' => $self->o('file_type')},
        'exclude_file_columns' => {},
        'require_file_attributes' => {},
        'exclude_file_attributes' => {},

        'tabix'	=> '/nfs/1000g-work/G1K/work/bin/tabix/tabix',
		'bgzip_exe'	=> '/nfs/1000g-work/G1K/work/bin/tabix/bgzip',
		'transpose_window_size' => 0,
		
		'depth_matrix_dir' 			=> '/nfs/1000g-work/G1K/work/zheng/annotate_p3_vcf/depth_matrix',  
		'supp_dp_mx_dir' 			=> '/nfs/1000g-work/G1K/work/zheng/annotate_p3_vcf/supp_depth_matrix',
		'sample_panel'				=> '/nfs/1000g-archive/vol1/ftp/release/20130502/integrated_call_samples.20130502.ALL.panel',
		'related_sample_list'		=> '/nfs/1000g-archive/vol1/ftp/release/20130502/20140625_related_individuals.txt',
		'simple_var_gl_file_dir'	=> '/nfs/1000g-work/G1K/scratch/zheng/add_gl_to_vcf',
		'complex_var_gl_file_dir'	=> '/nfs/1000g-archive/vol1/ftp/technical/working/20130723_phase3_wg/mvncall',

        'reference' => '/nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/bwa/grc37/human_g1k_v37.fa',        
        'fai' => $self->o('reference') . '.fai',
        
        final_output_dir => $self->o('root_output_dir'),
        #final_output_layout => '#final_output_dir#/#vcf_base_name#',
        name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method => 'basic',
        name_file_params => {
            new_dir => $self->o('final_output_dir'),
#           new_basename => '#vcf_base_name#',
#           add_datestamp => 1,
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
    	dir_label_params => ['vcf_base_name'], ### this is to decide the output dir structure in $self->output_dir
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '200Mb' => { 'LSF' => '-C0 -M200 -q '.$self->o('lsf_queue').' -R"select[mem>200] rusage[mem=200]"' },
            '500Mb' => { 'LSF' => '-C0 -M500 -q '.$self->o('lsf_queue').' -R"select[mem>500] rusage[mem=500]"' },
            '1Gb'   => { 'LSF' => '-C0 -M1000 -q '.$self->o('lsf_queue').' -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q '.$self->o('lsf_queue').' -R"select[mem>2000] rusage[mem=2000]"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q '.$self->o('lsf_queue').' -R"select[mem>4000] rusage[mem=4000]"' },
            '5Gb' => { 'LSF' => '-C0 -M5000 -q '.$self->o('lsf_queue').' -R"select[mem>5000] rusage[mem=5000]"' },
            '8Gb' => { 'LSF' => '-C0 -M8000 -q '.$self->o('lsf_queue').' -R"select[mem>8000] rusage[mem=8000]"' },
            '16Gb' => { 'LSF' => '-C0 -M16000 -q '.$self->o('lsf_queue').' -R"select[mem>16000] rusage[mem=16000]"' },
    };
}

sub pipeline_analyses {
    my ($self) = @_;
    my @analyses;
    push(@analyses, {
            -logic_name		=> 'get_seeds',
            -module			=> 'ReseqTrack::Hive::Process::SeedFactory',
            -meadow_type	=> 'LOCAL',
            -parameters		=> {
                seeding_module => $self->o('seeding_module'),
                seeding_options => $self->o('seeding_options'),
            },
            -flow_into => {
            	2 => { 	'generate_bed'  },
            },
      });     
  
     push(@analyses, {
          -logic_name    => 'generate_bed',
          -module        => 'ReseqTrack::Hive::Process::VcfToBed',
          -parameters    => {
              vcf => '#name#',
              bgzip => $self->o('bgzip_exe'),
              max_variants => $self->o('max_variants'),
          },
          -rc_name => '1Gb',
          -hive_capacity  =>  200,
          -flow_into => {
                '1' => { 'regions_factory_1' => { 	'vcf_bed' => '#bed#', 
                									'CHROM_start' => '#CHROM_start#', 
                									'CHROM_end' => '#CHROM_end#',
                									'input_vcf' => '#name#',
                 									'vcf_base_name' => '#expr( ($vcf =~ /([^\/]*)\.vcf(\.gz)?$/)[0] )expr#' }}
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
                 '2->A' =>  { 'annotate_vcf' => { 	'input_vcf' => '#input_vcf#',
                 									'region1' => '#vcf_base_name#.#SQ_start#.#bp_start#.#SQ_end#.#bp_end#',
                              						'SQ_start' => '#SQ_start#', 
                              						'bp_start' => '#bp_start#',
                              						'SQ_end' => '#SQ_end#', 
                              						'bp_end' => '#bp_end#',
                              						'fan_index' => '#fan_index#' },
                             },
                 'A->1' =>  [ 'merge_vcf' ],
            },
      }); 
	push(@analyses, {
            -logic_name    	=> 'annotate_vcf',
            -module        	=> 'ReseqTrack::Hive::Process::AnnotateVcf',
            -parameters		=> {
                tabix					=> $self->o('tabix'),
				depth_matrix_dir 		=> $self->o('depth_matrix_dir'),
				supp_dp_mx_dir			=> $self->o('supp_dp_mx_dir'),
				sample_panel 			=> $self->o('sample_panel'),
				related_sample_list		=> $self->o('related_sample_list'),
				simple_var_gl_file_dir	=> $self->o('simple_var_gl_file_dir'),
				complex_var_gl_file_dir	=> $self->o('complex_var_gl_file_dir'),
                output_param 			=> 'vcf',
            },
            -rc_name 		=> '2Gb',
          	-hive_capacity  =>  100,
          	-flow_into 		=> {
            		1 => [ ':////accu?vcf=[fan_index]', ':////accu?bp_start=[fan_index]', ':////accu?bp_end=[fan_index]' ],
          },
      });        
            
    push(@analyses, {
          -logic_name    => 'merge_vcf',
          -module        => 'ReseqTrack::Hive::Process::MergeVcf',
          -parameters    => {
              run_tabix		=> 1,
              bgzip => $self->o('bgzip_exe'),
              reseqtrack_options => {
                decode_file_id => 'vcf',
                #denestify => ['vcf','bp_start','bp_end'],
                #delete_param => ['vcf', 'vcf_bed'],   ### FIXME. add this line back after testing
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
              type => 'ANNOTATED_VCF',
              file => ['#vcf#'],
              name_file_module => $self->o('name_file_module'), 
              name_file_method => $self->o('name_file_method'),
              name_file_params => $self->o('name_file_params'),
            },
            -meadow_type => 'LOCAL', 
            -rc_name => '200Mb',
            -flow_into => {
              1 => ['mark_seed_complete'],
            }    
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

