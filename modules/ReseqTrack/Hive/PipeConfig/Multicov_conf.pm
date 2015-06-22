=head1 NAME

 ReseqTrack::Hive::PipeConfig::Multicov_conf

=head1 SYNOPSIS

  Pipeline must be seeded by the collection table of a ReseqTrack database (collection of genotype VCF files, one chromosome per collection)
  A multicov matrix .txt file will be created for each collection (chromosome), containing coverage info for all samples, and stored in the ReseqTrack database

  Here is an example pipeline configuration to load using reseqtrack/scripts/pipeline/load_pipeline_from_conf.pl

[Multicov]
table_name=file
config_options=-file_type VCF
config_options=-bam_type EXOME
config_module=ReseqTrack::Hive::PipeConfig::Multicov_conf
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


package ReseqTrack::Hive::PipeConfig::Multicov_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },

        'pipeline_name' => 'multicov',

        seeding_module => 'ReseqTrack::Hive::PipeSeed::BasePipeSeed',
        seeding_options => {
            output_columns => ['name', 'file_id'], 
            require_columns => $self->o('require_file_columns'),
            exclude_columns => $self->o('exclude_file_columns'),
            require_attributes => $self->o('require_file_attributes'),
            exclude_attributes => $self->o('exclude_file_attributes'),
          },
          
        'tabix'	=> '/nfs/1000g-work/G1K/work/bin/tabix/tabix',
		'BedTools_program' => '/nfs/1000g-work/G1K/work/bin/BEDTools/bin/multiBamCov',
		'paste_exe'	=> '/usr/bin/paste',
		
        'require_file_columns' => {'type' => $self->o('file_type')},
        'exclude_file_columns' => {},
        'require_file_attributes' => {},
        'exclude_file_attributes' => {},


        final_output_dir => $self->o('root_output_dir'),
        #final_output_layout => '#final_output_dir#/#vcf_base_name#',
        name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method => 'basic',
        name_file_params => {
            new_dir => $self->o('final_output_dir'),
            new_basename => '#vcf_base_name#',
            add_datestamp => 1,
            suffix => ['.matrix'],
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
        dir_label_params => ['vcf_base_name', 'sample'],  			### this is to decide the output dir structure in $self->output_dir
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
                2=> [ 'make_vcf_site_file'],     
            },
      });          
     push(@analyses, {
            -logic_name		=> 'make_vcf_site_file',
            -module			=> 'ReseqTrack::Hive::Process::GtVcf2SiteVcf', 
            -parameters		=> {
                output_param => 'site_vcf',
                gt_vcf => '#name#',  #### output_columns in seeding_option defines the seeds 
            },
            -flow_into => {
                1 => {'sample_factory' => {'site_vcf'=>'#site_vcf#', 'vcf_base_name' => '#expr( ($gt_vcf =~ /([^\/]*)\.vcf(\.gz)?$/)[0] )expr#'}},
                #1 => {'sample_factory' => {'site_vcf'=>'#site_vcf#', 'vcf_base_name' => '#expr( ($gt_vcf =~ /([^\/]*)\.biallelic_svm_snps_indelsAF0\.005_svs_indels_complex_STRs_svs\.sites\.vcf(\.gz)?$/)[0] )expr#'}},
            	### FIXME - use the first line for more generic uses
            },
      }); 
    push(@analyses, {
            -logic_name    => 'sample_factory',
            -module        => 'ReseqTrack::Hive::Process::SampleFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                tabix => $self->o('tabix'),
                #vcf => '#name#',  ### This is for situation when a full genotype VCF is provided as seed, rather than just the site files
                vcf => '/nfs/1000g-archive/vol1/ftp/technical/working/20130723_phase3_wg/shapeit2/ALL.chr22.phase3_shapeit2_integrated.20130502.snps_indels_svs.genotypes.vcf.gz', 
            	### FIXME: change to the upper line in other cases
            },
            -flow_into => {
                '2->A' => {'find_bam' => {'sample' => '#sample#', 'fan_index' => '#fan_index#', 'site_vcf'=>'#site_vcf#' }},
                'A->1' => ['merge_sample_depth'],
            },
      });    
      push(@analyses, {
            -logic_name    => 'find_bam',
            -module        => 'ReseqTrack::Hive::Process::FindBam',   
            -meadow_type => 'LOCAL',      
            -parameters    => {
                bam_type => $self->o('bam_type'),
                output_param => 'bam',
            },	
  			-flow_into => {
                1 => { 'calculate_depth' => {'bams'=>'#bam#', 'site_vcf'=>'#site_vcf#', 'fan_index' => '#fan_index#' } },
            },
      });
    push(@analyses, {
            -logic_name    => 'calculate_depth',
            -module        => 'ReseqTrack::Hive::Process::Multicov',         
            -parameters    => {
                program => $self->o('BedTools_program'),
                sample => '#sample#',
                stream_out => 0,
                output_param => 'sample_depth_file',
            },
            -hive_capacity  =>  200,
  			-flow_into => {
                1 => [':////accu?sample_depth_file=[fan_index]'],
            },
      });       
    push(@analyses, {
          -logic_name    => 'merge_sample_depth',
          -module        => 'ReseqTrack::Hive::Process::MergeByPaste',  
          -parameters    => {
              program => $self->o('paste_exe'),
              files => '#sample_depth_file#',
              site_vcf => '#site_vcf#',
              output_param => 'matrix',
              reseqtrack_options => {
                delete_param => ['sample_depth_file', 'site_vcf'],
              }
          },
          -rc_name => '16Gb',
          -analysis_capacity  =>  8,
          -flow_into	=> {
              1 => ['store_matrix'],
          }    
      });     
      push(@analyses, {
            -logic_name    => 'store_matrix',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -parameters    => {
              type => 'matrix',
              file => '#matrix#',
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

