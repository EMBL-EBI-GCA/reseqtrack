
=pod

=head1 NAME

ReseqTrack::Tools::BedTools_multicov

=head1 SYNOPSIS

This is a wrapper to run BedTools' multicov function to calculate read coverage.

It takes a VCF file (or BED, GFF file) that defines the genomic intervals of interest, it then calculate read depth from BAMs provided as input to the -bams option 
Can use multiple -bams tags to input multiple BAMs

The depth calculated from each BAM will be listed at the end of each line representing each genomic site in the VCF file, in the same order of the order of the input BAMs.

If -stream_out is used, the output will be written to STDOUT; otherwise, an output file with suffix .depth will be created 

example

my $multicov = $IR->new(
     -bams        =>"/path/to/bam1 /path/to/bam2" ,
     -bed        =>"/path/to/vcf/or/bed/file",
     -stream_out	=> 0,
     -program       =>"/path/to/multicov/executable",
     -working_dir     => '/path/to/dir/',
     -save_files_from_deletion => 1,
);

=cut

package ReseqTrack::Tools::BedTools_multicov;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);
use ReseqTrack::Tools::RunProgram;

use vars qw(@ISA);
@ISA = qw(ReseqTrack::Tools::RunProgram);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ($bed, $stream_out) = rearrange(
		[
			qw(
			  BED
			  STREAM_OUT
			  )
		],
		@args
	);
	
	$self->bed($bed);
	$self->stream_out($stream_out);
	return $self;
}

sub run_program {
	my $self = shift;
	
	#print "input bams are " ;
	#print join (" ", @{$self->input_files}) . "\n";
	#print "bed file is " . $self->bed . "\n";
	
	check_file_exists($_) foreach (@{$self->input_files});
	
	my @sh_input_names = ();
	if (@{$self->input_files} > 100) {
		my @sorted_input_files = sort {$a cmp $b} @{$self->input_files};
		my $short_name_hash = $self->get_short_input_names(2);
		foreach my $long_input_name ( @sorted_input_files ) {
			push @sh_input_names, $short_name_hash->{$long_input_name};
		}
	}	
	
	check_file_exists($self->bed);
	check_executable($self->program);

	my @cmd_words = ($self->program);
	
	if (@{$self->input_files} > 100) {
		push(@cmd_words, '-bams', join(" ", @sh_input_names));
	}
	else {
		push(@cmd_words, '-bams', join(" ", @{$self->input_files}));
	}
	push(@cmd_words, '-bed', $self->bed);
	#push(@cmd_words, '-r', $self->options('r')) if ($self->options('r')); ## the -r function is not supported
	
	unless ($self->stream_out) {
		my $out_vcf = $self->working_dir . "/" . basename($self->bed) . "." . $self->job_name . ".$$.depth";  
		push (@cmd_words, ">", $out_vcf);
		$self->output_files($out_vcf);
	}
	
	my $cmd = join(' ', @cmd_words);

	$self->execute_command_line ($cmd);

	return;

}	

sub bed {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{bed} = $arg;
  }
  return $self->{bed};
}

sub stream_out {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{stream_out} = $arg;
  }
  return $self->{stream_out};
}