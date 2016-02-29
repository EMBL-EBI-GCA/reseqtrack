package ReseqTrack::Hive::PipeConfig::WGBS_conf;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub default_options {
    my ($self) = @_;

    return {
        %{ $self->SUPER::default_options() },
        'fastqc_exe' => '/hps/cstor01/nobackup/faang/ernesto/bin/FastQC/fastqc',
	'fastqscreen_exe' => '/hps/cstor01/nobackup/faang/ernesto/bin/FastQScreen_v0.3.1/fastq_screen',
	'fastqscreen_conf' => '/hps/cstor01/nobackup/faang/ernesto/reference/fastq_screen_databases/fastq_screen.conf',
	'split_exe' => $self->o('ENV', 'RESEQTRACK').'/c_code/split/split',
	'chunk_max_reads' => 5000000,
	 'bwa_exe' => '/nfs/1000g-work/G1K/work/bin/bwa/bwa',
        'samtools_exe' => '/nfs/1000g-work/G1K/work/bin/samtools/samtools',
    };
}


sub pipeline_analyses {
    my ($self) = @_;

    return [
        {
            -logic_name => 'run_fastqc',
            -module     => 'ReseqTrack::Hive::Process::RunFastqc',
            -parameters => {
                'fastq'      => '#filename#',
                'output_dir' => '#output_dir#',
                'program_file' => $self->o('fastqc_exe')
            },
            -flow_into => {
                1 => {'run_fastqscreen'=> {'fastq'=> '#filename#', 'output_dir' => '#output_dir#'}},
            },
        },
        {
            -logic_name => 'run_fastqscreen',
            -module     => 'ReseqTrack::Hive::Process::RunFastQScreen',
            -parameters => {
                'output_dir' => '#output_dir#',
                'program_file' => $self->o('fastqscreen_exe'),
                'conf_file' => $self->o('fastqscreen_conf')
            },
	    -flow_into => {
		 1 =>  {'split_fastq'=> {'fastq'=> '#fastq#','output_dir'=> '#output_dir#'}}
	    },
        },
	{
	    -logic_name    => 'split_fastq',
            -module        => 'ReseqTrack::Hive::Process::SplitFastq',
            -parameters    => {
		program_file => $self->o('split_exe'),
	        max_reads => $self->o('chunk_max_reads'),
		'run_source_id' => 'input1000'
            },
	}
    ];
}

1;

