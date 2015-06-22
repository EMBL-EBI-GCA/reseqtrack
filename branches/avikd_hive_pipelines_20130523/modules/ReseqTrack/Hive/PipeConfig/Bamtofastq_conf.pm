=head1 NAME

 ReseqTrack::Hive::PipeConfig::Bamtofastq_conf

=head1 SYNOPSIS
  
  Options you MUST specify on the command line:
  -study_id, refers to a study_id in your run_meta_info table. Can be specified multiple times.
  -password, for accessing the hive database
  -reseqtrack_db_name, (or -reseqtrack_db -db_name=??) your reseqtrack database 

 Connection to the hive database:
  -pipeline_db -host=???, (default )
  -pipeline_db -port=???, (default 4194)
  -pipeline_db -user=???, must have write access (default )
  -dipeline_db -dbname=???, (default is a mixture of your unix user name + the pipeline name)

 Connection to the reseqtrack database:
  -reseqtrack_db -host=???, 
  -reseqtrack_db -user=???, 
  -reseqtrack_db -port=???, 
  -reseqtrack_db -pass=???, 

  -root_output_dir, (default is your current directory)
  -type_bam, type of BAM files to look for in the reseqtrack database, default EGA_BAM
  -type_fastq, output type for fastq files in the reseqtrack database, default EGA_FASTQ

=cut

package ReseqTrack::Hive::PipeConfig::Bamtofastq_conf;

use strict;
use warnings;

use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },
        'pipeline_name' => 'bamtofastq',
        'type_bam'    => 'EGA_BAM',
        'type_fastq'    => 'EGA_FASTQ',
        'picard_dir' => '/nfs/ega/private/ega/work/ext/Tools/picard-tools-1.88',
        'study_id' => [],
       };
}

sub pipeline_create_commands {
    my ($self) = @_;

    return [
        @{$self->SUPER::pipeline_create_commands},
    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '500Mb' => { 'LSF' => '-C0 -M500 -q analysis -R"select[mem>500] rusage[mem=500]"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q analysis -R"select[mem>4000] rusage[mem=4000]"' },
            '6Gb' => { 'LSF' => '-C0 -M6000 -q analysis -R"select[mem>6000] rusage[mem=6000]"' },
        };
}

sub pipeline_analyses {
    my ($self) = @_;

    my @analyses;
    push(@analyses, {
            -logic_name    => 'studies_factory',
            -module        => 'ReseqTrack::Hive::Process::JobFactory',
            -meadow_type => 'LOCAL',
            -input_ids => [{study_id => $self->o('study_id')}],
            -parameters    => {
                factory_value => '#study_id#',
                temp_param_sub => { 2 => [['study_id','factory_value']]}, # temporary hack pending updates to hive code
            },
            -flow_into => {
                2 => [ 'runs_factory' ],
            },
      });
       push(@analyses, {
            -logic_name    => 'runs_factory',
            -module        => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                type_branch => 'run',
            },
            -flow_into => {
		2 => [ 'library_layout_factory' ],
          },
      });
 	push(@analyses, {
            -logic_name    => 'library_layout_factory',
            -module        => 'ReseqTrack::Hive::Process::RunMetaInfoFactory',
            -meadow_type => 'LOCAL',
            -parameters    => {
                type_branch => 'library_layout',
            },
            -flow_into => {
		'2->A' => [ 'find_source_bams' ],
                'A->1' => ['final_step'],
           },
      });
	push(@analyses, {
            -logic_name    => 'find_source_bams',
            -module        => 'ReseqTrack::Hive::Process::ImportCollection',
            -meadow_type => 'LOCAL',
            -parameters    => {
                collection_type => $self->o('type_bam'),
                collection_name => '#run_id#',
                output_param => 'bam',
            },
            -flow_into => {
                1 => [ 'revert_bam', ':////accu?bam=[]' ],
            },
      });
	push(@analyses, {
            -logic_name => 'revert_bam',
            -module        => 'ReseqTrack::Hive::Process::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                command => 'revert_sam',
                create_index => 1,
                options     => {validation_stringency => 'SILENT', remove_duplicate_info => 1, remove_alignment_info => 1,},
                jvm_args => '-Xmx4g',
                
            },
            -rc_name => '4Gb',
            -hive_capacity  =>  200,
  	    -analysis_capacity => 100,
            -flow_into => {
                1 => [ 'bam_to_fastq'],
            },
      });
	push(@analyses, {
            -logic_name => 'bam_to_fastq',
            -module        => 'ReseqTrack::Hive::Process::RunPicard',
            -parameters => {
                picard_dir => $self->o('picard_dir'),
                command => 'sam_to_fastq',
                jvm_args => '-Xmx4g',
                run_id => '#run_id#',
                output_param => 'fastq',
                library_layout => '#library_layout#',
            },
            -rc_name => '4Gb',
            -hive_capacity  =>  200,            
  	    -analysis_capacity => 100,
	    -flow_into => {
                1 => ['load_fastq'],
            },  
      });
	push(@analyses, {
            -logic_name    => 'load_fastq',
            -module        => 'ReseqTrack::Hive::Process::LoadFile',
            -parameters    => {
                type => $self->o('type_fastq'),
                file => '#fastq#',
                collection_name => '#run_id#',
		collect => 1,
 	        delete_param => ['bam','bai'],
            },
            -rc_name => '500Mb',
  	    -analysis_capacity => 30,
      });
	push(@analyses, {
            -logic_name    =>'final_step',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
  	    -meadow_type => 'LOCAL',
           -parameters  => {
                        "cmd" => "echo done #run_id# ",
              },
          });

    return \@analyses;
}

1;
