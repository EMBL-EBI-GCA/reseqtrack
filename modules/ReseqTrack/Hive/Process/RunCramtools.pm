
package ReseqTrack::Hive::Process::RunCramtools;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use ReseqTrack::Tools::ConvertBam;
use File::Basename;

sub param_defaults {
  return {
    cramtools_options => {},
    program_file => undef,
    reference => undef,
  };
}


sub run {
    #throw("fixing cramtools3...."); 
    my $self = shift @_;
    $self->param_required('bam');
    my $bams = $self->param_as_array('bam');
	#$self->param_required('reference');
 	
    my $options = $self->param('cramtools_options');

	my $output_file_name = $self->output_dir . "/" . basename($bams->[0]) . ".cram";
	my $cram = ReseqTrack::Tools::ConvertBam->new (
		-input_files 		=> $bams,
		-program			=> $self->param('program_file'),  
		-cramtools_jar_file => $self->param('cramtools_jar_file'),
		-java_exe			=> $self->param('java_exe'),
		-jvm_args			=> $self->param('jvm_args'),
		-output_format 		=> 'cram',
		-options 			=> $options,
		-reference_fasta 	=> $self->param('reference'),
		-working_dir		=> $self->output_dir,
		-output_file 		=> $output_file_name,
	);

    $self->run_program($cram);
    my $output_files = $cram->output_files;
	
	if (basename($output_files->[0]) =~ /\.cram$/i) {
		$self->output_param('cram'  => $output_files);
	}
	elsif (basename($output_files->[0]) =~ /\.crai$/i) {
		$self->output_param('crai'  => $output_files);
	}
		
}


1;

