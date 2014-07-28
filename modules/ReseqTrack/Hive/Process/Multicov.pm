package ReseqTrack::Hive::Process::Multicov;

use strict;
use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use File::Path qw(mkpath);
use File::Basename qw(basename);
use ReseqTrack::Tools::BedTools_multicov;
use IPC::System::Simple qw(system);

sub param_defaults {
  return {
    bams => undef,
    site_vcf => undef,
    sample => undef,
    program => undef,
    stream_out => undef,
    options => undef,
    save_files_from_deletion => 0,
  };
}

sub run {
    my $self = shift @_;

    $self->param_required('bams');
    my $bams =  $self->param_as_array('bams');   
    
    my $site_vcf = $self->param_required('site_vcf');
	my $sample = $self->param('sample');
	
    my $program = $self->param('program');
	my $stream_out = $self->param('stream_out');
	my $options = $self->param('options');
	my $save_files_from_deletion = $self->param('save_files_from_deletion');

 	my $run_multicov = ReseqTrack::Tools::BedTools_multicov->new(
 		-program					=> $program,
		-input_files				=> $bams,
		-bed						=> $site_vcf,
		-stream_out					=> $stream_out,
		-working_dir				=> $self->output_dir,
		-options					=> $options,
		-job_name					=> $sample,
		-save_files_from_deletion	=> $save_files_from_deletion,
	);

	$self->run_program($run_multicov);
    my $output_files = $run_multicov->output_files;
    
    mkpath($self->output_dir) unless (-d $self->output_dir);
    
    my $tmp_out = $output_files->[0] . ".tmp";
    my $cut_cmd = "cut -f1-8 --complement " . $output_files->[0] . " > " . $tmp_out;  ## This is to get the last columns of the vcf file which contain depth
    system($cut_cmd);
    my $exit = $?>>8;
	throw("cut VCF $cut_cmd failed\n") if ($exit >=1);
    
    my $header = $self->output_dir . "/" .   basename($site_vcf) . "." . $sample . ".$$.header";
	open(OUT, ">", $header) || throw("Cannot open out file $header");
	print OUT "$sample\n";
	close OUT;
	
	my $depth_only_file = $tmp_out . ".depth_only";
	system("cat $header $tmp_out > $depth_only_file");
	my $exit2 = $?>>8;
	throw("cat files failed\n") if ($exit2 >=1);

    $self->output_param('sample_depth_file', $depth_only_file);
    
    unlink($tmp_out);
    unlink($header);
    unlink($output_files->[0]);
}

1;