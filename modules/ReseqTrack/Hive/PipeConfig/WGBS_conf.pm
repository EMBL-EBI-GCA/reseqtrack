package ReseqTrack::Hive::PipeConfig::WGBS_conf;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },
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
	'RGSM'                => 'test_s',
        'RGPU'                => 'test_r'
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

    return [
	{
            -logic_name => 'find_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                'inputcmd'      => 'find #directory# -type f', 
		'column_names'  => ['fastq']
            },
	    -input_ids => [
		 {
		     'directory' => '/hps/cstor01/nobackup/faang/ernesto/wgbs_pipeline/data/',
		     'o_dir' => '/hps/cstor01/nobackup/faang/ernesto/wgbs_pipeline/bismark_out/',
		     'run_id' => 'test_data'
		 }
	    ],
            -flow_into => {
		2 => ['split_fastq'],
#		2 => ['run_fastqc','run_fastqscreen', 'split_fastq'],
            },
        },
	{
           -logic_name => 'run_fastqc',
           -module     => 'ReseqTrack::Hive::Process::RunFastqc',
           -parameters => {
               'fastq'      => '#directory# #fastq#',
               'output_dir' => '#o_dir#',
               'program_file' => $self->o('fastqc_exe')
            },
	    -analysis_capacity => 4,
	    -rc_name => '200Mb'
        },
	{
            -logic_name => 'run_fastqscreen',
            -module     => 'ReseqTrack::Hive::Process::RunFastQScreen',
            -parameters => {
                'fastq'      => '#fastq#',
                'output_dir' => '#o_dir#',
                'program_file' => $self->o('fastqscreen_exe'),
                'conf_file' => $self->o('fastqscreen_conf')
             },
	    -analysis_capacity => 4,
	    -rc_name => '1Gb'
        },
	{
           -logic_name    => 'split_fastq',
           -module        => 'ReseqTrack::Hive::Process::SplitFastq',
           -parameters    => {
               'fastq'      => '#fastq#',
               'output_dir' => '#o_dir#',
               program_file => $self->o('split_exe'),
               max_reads => $self->o('chunk_max_reads'),
               'run_source_id' => '#run_id#',
            },
	   -flow_into => {
	      '2->A' => {'bismark_mapper' => {'chunk_label' => '#run_source_id#.#chunk#', 'fastq' => '#fastq#'}},
              'A->1' => ['merge_chunks'],
	   },
	    -analysis_capacity => 4,
	    -rc_name => '200Mb',
        },
	{
	    -logic_name => 'bismark_mapper',
            -module        => 'ReseqTrack::Hive::Process::RunBismark',
            -parameters    => {
                program_file => $self->o('bismark_exe'),
                samtools => $self->o('samtools_exe'),
                reference => $self->o('reference'),
		command => 'aln',
	        RGSM => $self->o('RGSM'),
                RGPU => $self->o('RGPU'),
		'output_dir' => '#o_dir#',
		run_id => '#run_id#'
            },
	    -flow_into => {
		1 => ['sort_chunks']
	    },
	     -analysis_capacity => 20,
	     -rc_name => '10Gb'
	},
	{
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
       },
       {
         -logic_name => 'merge_chunks',
         -module        => 'ReseqTrack::Hive::Process::RunBioBamBam',
         -parameters => {
             biobambam_dir => $self->o('biobambam_dir'),
             command => 'bammerge',
             create_index => 1,
	     'chunk_label' => '#run_id#',
	     'output_dir' => '#o_dir#'
	 },
	 -flow_into => {
              1 => [ 'mark_duplicates'],
	 },
        },
	{
	    -logic_name => 'mark_duplicates',
            -module        => 'ReseqTrack::Hive::Process::RunBioBamBam',
	    -parameters => {
		biobambam_dir => $self->o('biobambam_dir'),
                command => 'bammarkduplicates',
		'output_dir' => '#o_dir#',
                create_index => 1,
		'chunk_label' => '#run_id#'
	    },
	    -flow_into => {
		1 => ['bismark_methcall'],
	    },
        },
	{
	     -logic_name => 'bismark_methcall',
             -module        => 'ReseqTrack::Hive::Process::RunBismark',
	     -parameters    => {
                program_file => $self->o('bismark_methcall_exe'),
                command => 'methext',
		'LIBRARY_LAYOUT' => 'SINGLE',
                'output_dir' => '#o_dir#'
	     },
	}
    ];
}

1;

