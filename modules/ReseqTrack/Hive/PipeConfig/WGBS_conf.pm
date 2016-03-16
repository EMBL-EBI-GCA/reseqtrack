package ReseqTrack::Hive::PipeConfig::WGBS_conf;

use strict;
use warnings;
use base ('ReseqTrack::Hive::PipeConfig::ReseqTrackGeneric_conf');                                                               

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },
	    seeding_module => 'ReseqTrack::Hive::PipeSeed::BasePipeSeed',
            seeding_options => {
		collection_type                => 'WGBS_FASTQ',
		output_columns                 => $self->o('run_columns'),
		output_attributes              => $self->o('run_attributes'),
		require_columns                => $self->o('require_run_columns'),
		exclude_columns                => $self->o('exclude_run_columns'),
		require_attributes             => $self->o('require_run_attributes'),
		exclude_attributes             => $self->o('exclude_run_attributes'),
		path_names_array               => $self->o('path_names_array'),
	    },

	 regexs               => undef,
	 path_names_array     => [ 'sample_desc_1', 'sample_desc_2', 'sample_desc_3', 'library_strategy', 'center_name' ],
 	 type_fastq           => 'WGBS_FASTQ',

	'biobambam_dir' => '/nfs/production/reseq-info/work/bin/biobambam2-2.0.10/bin/',
        'fastqc_exe' => '/hps/cstor01/nobackup/faang/ernesto/bin/FastQC/fastqc',
	'fastqscreen_exe' => '/hps/cstor01/nobackup/faang/ernesto/bin/FastQScreen_v0.3.1/fastq_screen',
	'fastqscreen_conf' => '/hps/cstor01/nobackup/faang/ernesto/reference/fastq_screen_databases/fastq_screen.conf',
	'lsf_queue' => 'production-rh6',
	'split_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/split/split',
	'chunk_max_reads' => 5000000,
	 'bismark_exe' => '/hps/cstor01/nobackup/faang/ernesto/bin/bismark_v0.15.0/bismark',
	'bismark_methcall_exe' => '/hps/cstor01/nobackup/faang/ernesto/bin/bismark_v0.15.0/bismark_methylation_extractor',
        'samtools_exe' => '/hps/cstor01/nobackup/faang/ernesto/bin/samtools-1.3/samtools',
	'reference' => '/hps/cstor01/nobackup/faang/ernesto/reference/bismark',

        'RGSM'                => '#sample_source_id#',
        'RGPU'                => '#run_source_id#',

	'bam_type'            => 'WGBS_RUN_BAM',
        'collection_name'     => '#run_source_id#',
        'build_collection'    => 1,

        'run_attributes'        => [],
        'run_columns'           => ['run_id', 'run_source_id', 'center_name', 'run_alias'],

	require_run_attributes        => {},
	exclude_run_attributes        => {},
	require_run_columns           => { status => ['public', 'private'], },
	exclude_run_columns           => {},

        'dir_label_params_list'       => ["sample_source_id", "experiment_source_id", "run_source_id", "chunk_label"],

        final_output_dir              => "/hps/cstor01/nobackup/faang/ernesto/wgbs_pipeline",
        final_output_layout           => '#sample_desc_1#/#sample_desc_2#/#sample_desc_3#/#library_strategy#/#center_name#',
        name_file_module              => 'ReseqTrack::Hive::NameFile::BaseNameFile',
        name_file_method              => 'basic',
	name_file_params              => {
                                          new_dir       => '#final_output_dir#/#final_output_layout#',
                                          new_basename  => '#run_source_id#.bwa.GRCh38',
                                          add_datestamp => 1,
                                          suffix         => '.bam',
                                         },
    };
}

sub resource_classes {
    my ($self) = @_;
    return {
            %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
            '200Mb' => { 'LSF' => '-C0 -M200 -q '.$self->o('lsf_queue').' -R"select[mem>200] rusage[mem=200]"' },
            '500Mb' => { 'LSF' => '-C0 -M500 -q '.$self->o('lsf_queue').' -R"select[mem>500] rusage[mem=500]"' },
            '800Mb' => { 'LSF' => '-C0 -M800 -q '.$self->o('lsf_queue').' -R"select[mem>800] rusage[mem=800]"' },
            '1Gb'   => { 'LSF' => '-C0 -M1000 -q '.$self->o('lsf_queue').' -R"select[mem>1000] rusage[mem=1000]"' },
            '2Gb' => { 'LSF' => '-C0 -M2000 -q '.$self->o('lsf_queue').' -R"select[mem>2000] rusage[mem=2000]"' },
            '3Gb' => { 'LSF' => '-C0 -M3000 -q '.$self->o('lsf_queue').' -R"select[mem>3000] rusage[mem=3000]"' },
            '4Gb' => { 'LSF' => '-C0 -M4000 -q '.$self->o('lsf_queue').' -R"select[mem>4000] rusage[mem=4000]"' },
            '5Gb' => { 'LSF' => '-C0 -M5000 -q '.$self->o('lsf_queue').' -R"select[mem>5000] rusage[mem=5000]"' },
            '6Gb' => { 'LSF' => '-C0 -M6000 -q '.$self->o('lsf_queue').' -R"select[mem>6000] rusage[mem=6000]"' },
            '8Gb' => { 'LSF' => '-C0 -M8000 -q '.$self->o('lsf_queue').' -R"select[mem>8000] rusage[mem=8000]"' },
	    '10Gb' => { 'LSF' => '-C0 -M10000 -q '.$self->o('lsf_queue').' -R"select[mem>10000] rusage[mem=10000]"' },
    };
}


sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
        %{$self->SUPER::pipeline_wide_parameters},

        dir_label_params => ["run_source_id", "chunk_label"],
    };
}

sub hive_meta_table {
    my ($self) = @_;
    return {
        %{$self->SUPER::hive_meta_table},       # here we inherit anything from the base class
        'hive_use_param_stack'  => 1,           # switch on the implicit parameter propagation mechanism
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
           2 => [ 'find_source_fastqs' ],
         },
     });
    
    push(@analyses, {
        -logic_name    => 'find_source_fastqs',
        -module        => 'ReseqTrack::Hive::Process::ImportCollection',
        -meadow_type => 'LOCAL',
        -parameters    => {
            collection_type => $self->o('type_fastq'),
            collection_name => '#run_source_id#',
            output_param => 'fastq',
            reseqtrack_options => {
              flows_do_count_param => 'fastq',
              flows_do_count => {
                 1 => '1+',
                 2 => '0',
	      },
            }
	},
	-flow_into => {
                '1' => [ 'split_fastq' ],
                '2' => [ 'mark_seed_futile' ],
	},
     });

    push(@analyses, {
       -logic_name    => 'split_fastq',
       -module        => 'ReseqTrack::Hive::Process::SplitFastq',
       -parameters    => {
           program_file => $self->o('split_exe'),
           max_reads    => $self->o('chunk_max_reads'),
           regexs       => $self->o('regexs'),
       },
       -rc_name => '200Mb',
       -analysis_capacity  =>  4,
       -hive_capacity  =>  200,
       -flow_into => {
	   '2->A' => {'bismark_mapper' => {'chunk_label' => '#run_source_id#.#chunk#', 'fastq' => '#fastq#'}},
           'A->1' => ['merge_chunks'],
       }
    });

    push(@analyses, {
       -logic_name => 'bismark_mapper',
       -module        => 'ReseqTrack::Hive::Process::RunBismark',
       -parameters    => {
           program_file => $self->o('bismark_exe'),
           samtools => $self->o('samtools_exe'),
           reference => $self->o('reference'),
           command => 'aln',
           RGSM => $self->o('RGSM'),
           RGPU => $self->o('RGPU'),
          'output_dir' => '#o_dir#'
       },
       -analysis_capacity => 20,
       -rc_name => '10Gb',
       -flow_into => {
           1 => ['sort_chunks'],
       },
    });

    push(@analyses, {
       -logic_name => 'sort_chunks',
       -module        => 'ReseqTrack::Hive::Process::RunBioBamBam',
       -parameters => {
            biobambam_dir => $self->o('biobambam_dir'),
            command => "bamsort",
            'output_dir' => '#o_dir#',
             options => { fixmates => 1 },
             create_index => 1
       },
       -analysis_capacity  =>  20,
       -flow_into => {
	   1 => [ ':////accu?bam=[]', ':////accu?bai=[]']
	},
     });

    push(@analyses, {
       -logic_name => 'merge_chunks',
       -module        => 'ReseqTrack::Hive::Process::RunBioBamBam',
       -parameters => {
          biobambam_dir => $self->o('biobambam_dir'),
          command => 'bammerge',
          create_index => 1,
          'chunk_label' => '#run_source_id#',
          'output_dir' => '#o_dir#'
       },
       -flow_into => {
          1 => [ 'store_bam']
       },
    });

    push(@analyses, {
        -logic_name    => 'store_bam',
        -module        => 'ReseqTrack::Hive::Process::LoadFile',
        -parameters    => {
            type => $self->o('bam_type'),
            file => '#bam#',
            name_file_module    => $self->o('name_file_module'),
            name_file_method    => $self->o('name_file_method'),
            name_file_params    => $self->o('name_file_params'),
            final_output_dir    => $self->o('final_output_dir'),
            final_output_layout => $self->o('final_output_layout'),
            collection_name     => $self->o('collection_name'),
            collect             => $self->o('build_collection'),
        },
        -rc_name => '200Mb',
        -hive_capacity  =>  200,
    });

    push(@analyses, { 
        -logic_name => 'bismark_methcall',
        -module        => 'ReseqTrack::Hive::Process::RunBismark',
	-parameters    => {
	    program_file => $self->o('bismark_methcall_exe'),
	    command => 'methext',
	    'output_dir' => '#o_dir#'
	},
    }); 

    push(@analyses, {
       -logic_name    => 'mark_seed_futile',
       -module        => 'ReseqTrack::Hive::Process::UpdateSeed',
       -parameters    => {
          is_futile  => 1,
       },
       -meadow_type => 'LOCAL',
   });
 
    return \@analyses;
}

1;

