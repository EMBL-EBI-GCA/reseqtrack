package ReseqTrack::Tools::QC::FastQC;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Copy qw(move);
use File::Basename qw(fileparse);

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

sub run_program {
	my ( $self, @args ) = @_;
	
	my $program = $self->program;
	my $temp_dir = $self->get_temp_dir;
	my @fastq_files = @{$self->input_files};
	
	my $cmd  = "$program -q -o $temp_dir ".join(' ',@fastq_files);
	
	$self->execute_command_line($cmd); 
        
        foreach my $fastq (@fastq_files) {
          print "Here1: $fastq\n";
          my $fastq_base_name = fileparse($fastq, qr/\.f(?:ast)?q(?:\.gz)?/);
          my $fastq_temp_dir = $temp_dir .'/'.$fastq_base_name .'_fastqc';
          my ($output_dir, $rename_files);

          if ($self->keep_file_output_dir) {
            $output_dir = $self->working_dir . '/'. $fastq_base_name . '_fastqc';
            $rename_files = 0;
          }
          else {
            $output_dir = $self->working_dir;
            $rename_files = 1;
          }
          $self->change_dir($output_dir); # to make sure directory exists

          if ($self->keep_zip) {
            my $from = "$fastq_temp_dir.zip";
            my $to = $output_dir.'/'.$fastq_base_name.'_fastqc.zip';
            $self->output_files($to);
            print "moving $from $to\n";
            move($from, $to) or throw("could not move $from to $to $!");
          }

          if ($self->keep_summary) {
            my $from = "$fastq_temp_dir/summary.txt";
            my $to = $rename_files ? $output_dir.'/'.$fastq_base_name.'_fastqc_summary.txt'
                    : "$output_dir/summary.txt";
            $self->output_files($to);
            print "moving $from $to\n";
            move($from, $to) or throw("could not move $from to $to $!");
          }

          if ($self->keep_text) {
            my $from = "$fastq_temp_dir/fastqc_data.txt";
            my $to = $rename_files ? $output_dir.'/'.$fastq_base_name.'_fastqc_report.txt'
                    : "$output_dir/fastqc_data.txt";
            $self->output_files($to);
            print "moving $from $to\n";
            move($from, $to) or throw("could not move $from to $to $!");
          }

          if ($self->keep_html) {
            foreach my $file_name ('fastqc_report.html', 'Icons', 'Images') {
              my $from = "$fastq_temp_dir/$file_name";
              my $to = $rename_files ? $output_dir.'/'.$fastq_base_name.'_'.$file_name
                      : "$output_dir/$file_name";
              $self->output_files($to);
              print "moving $from $to\n";
              move($from, $to) or throw("could not move $from to $to $!");
            }
          }
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
