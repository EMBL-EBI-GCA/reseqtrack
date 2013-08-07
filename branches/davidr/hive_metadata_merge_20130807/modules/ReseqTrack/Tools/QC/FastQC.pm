package ReseqTrack::Tools::QC::FastQC;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use base qw(ReseqTrack::Tools::RunProgram);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $keep_html, $keep_zip, $keep_text, $keep_summary, $keep_output_dir )	= rearrange( [
	    qw( KEEP_HTML KEEP_ZIP KEEP_TEXT KEEP_SUMMARY KEEP_OUTPUT_DIR )
		], @args);

	$self->keep_html($keep_html || 0);
	$self->keep_zip($keep_zip || 0);
	$self->keep_text($keep_text || 0);
	$self->keep_summary($keep_summary || 0);
	$self->keep_file_output_dir($keep_output_dir || 0);

	return $self;	
}

sub run {
	my ( $self, @args ) = @_;
	
	my $program = $self->program;
	my $working_dir = $self->working_dir;
	my @fastq_files = @{$self->input_files};
	
	my $cmd  = "$program -q -o $working_dir ".join(' ',@fastq_files);
	
	$self->execute_command_line($cmd); 
	
	for my $fastq_file (@fastq_files){
		$self->handle_output($fastq_file);
	}
}

sub handle_output {
	my ($self, $fastq_file) = @_;
	
	$self->keep_or_delete($self->keep_html, $self->html_report_paths($fastq_file));
	$self->keep_or_delete($self->keep_text, $self->report_text_path($fastq_file));
	$self->keep_or_delete($self->keep_summary, $self->summary_text_path($fastq_file));
	$self->keep_or_delete($self->keep_zip, $self->zipped_output_path($fastq_file));
	$self->keep_or_delete($self->keep_file_output_dir, $self->file_output_path($fastq_file));
}

sub output_base_name {
	my ($self, $fastq_file) = @_;
	my @chunks =  split /\//, $fastq_file;
	my $base_name = $chunks[-1];
	
	$base_name =~ s/\.fastq.gz/_fastqc/;
	
	my $working_dir = $self->working_dir;
	my $output_dir = "$working_dir/$base_name";
	
	return ($base_name,$output_dir);
}

sub file_output_path {
	my ($self, $fastq_file) = @_;
	$self->output_base_name($fastq_file);
	my ($base_name,$output_dir) = $self->output_base_name($fastq_file);	
	return $output_dir;
}

sub summary_text_path {
	my ($self, $fastq_file) = @_;
	$self->output_base_name($fastq_file);
	my ($base_name,$output_dir) = $self->output_base_name($fastq_file);	
	return "$output_dir/summary.txt";
}

sub report_text_path {
	my ($self, $fastq_file) = @_;
	$self->output_base_name($fastq_file);
	my ($base_name,$output_dir) = $self->output_base_name($fastq_file);	
	return "$output_dir/fastqc_data.txt";
}

sub html_report_paths {
	my ($self, $fastq_file) = @_;
	$self->output_base_name($fastq_file);
	my ($base_name,$output_dir) = $self->output_base_name($fastq_file);
	return ("$output_dir/fastqc_report.html","$output_dir/Icons","$output_dir/Images");
}

sub zipped_output_path {
	my ($self, $fastq_file) = @_;
	$self->output_base_name($fastq_file);
	my ($base_name,$output_dir) = $self->output_base_name($fastq_file);	
	my $working_dir = $self->working_dir;
	return "$working_dir/$base_name.zip";
}

sub keep_or_delete {
	my ($self,$keep,@files) = @_;
	
	if ($keep){
		$self->output_files(\@files);
	}
	else {
		$self->created_files(\@files);
	}
	
}

sub keep_file_output_dir {
	my ( $self, $arg ) = @_;
	  if ($arg) {
	    $self->{'keep_file_output_dir'} = $arg;
	  }
	  return $self->{'keep_file_output_dir'};
}

sub keep_html {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'keep_html'} = $arg;
  }
  return $self->{'keep_html'};
}

sub keep_zip {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'keep_zip'} = $arg;
  }
  return $self->{'keep_zip'};
}

sub keep_text {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'keep_text'} = $arg;
  }
  return $self->{'keep_text'};
}

sub keep_summary {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'keep_summary'} = $arg;
  }
  return $self->{'keep_summary'};
}

1;