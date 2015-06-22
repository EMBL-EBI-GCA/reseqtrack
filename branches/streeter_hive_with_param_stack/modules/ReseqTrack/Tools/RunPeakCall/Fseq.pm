package ReseqTrack::Tools::RunPeakCall::Fseq;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunProgram;

use base qw(ReseqTrack::Tools::RunPeakCall);

=pod

=head1 NAME

ReseqTrack::Tools::RunPeakCall::Fseq

=head1 SYNOPSIS

class for running F-seq. Child class of ReseqTrack::Tools::RunPeakCall
http://fureylab.web.unc.edu/software/fseq/

Set fragment_size to 0 for DNAse data

Takes either a bed file or a bam file.
BAM files will be converted to bed with bamToBed
If the BAM file has duplicates marked, the module will remove them prior to conversion if strip_duplicates is set to a value that evaluates to true

=head1 Example

my $fseq = ReseqTrack::Tools::RunPeakCall::Fseq(
					-input_files => '/path/to/file.bed'
					-working_dir => $output_dir,
					-options => \%options,
					-job_name => $name,
					-strip_duplicates => 1,
                      );
$fseq->run;
my $output_file_list = $fseq->output_files;

=cut
sub DEFAULT_OPTIONS { return {
		background_dir => undef, #	-b <background dir>     background directory (default=none)
		read_count => undef, #	 -c <arg>                genomic count of sequence reads (defualt = calculated)
		fragment_size => undef, # -f <arg>                fragment size (default=estimated from data)
		feature_length => undef, # -l <arg>                feature length (default=600)
		output_format => 'bed',  #-of <wig | bed | npf>   output format (default wig)
		ploidy_dir => undef, #	 -p <ploidy dir>         ploidy/input directory (default=none)
		track_step => undef, # -s <arg>                wiggle track step (default=1)
		threshold => undef, #	 -t <arg>                threshold (standard deviations) (default=4.0)
		verbose => 0, #	 -v                      verbose output
		wg_threshold => undef, #  -wg <arg>               wg threshold set (defualt = calculated)
	};
}

sub new {
	my ( $class, @args ) = @_;
	
	my $self = $class->SUPER::new(@args);

	#setting defaults
	if (!$self->program) {
	   if ($ENV{fseq}) {
	     $self->program($ENV{fseq} . '/fseq');
	   }
	   else {
	     $self->program('fseq');
	   }
	 }

	return $self;
}

sub can_read_bam {
	return 0;
}

sub can_use_control {
	return 0;
}
sub control_required {
	return 0;
}

sub run_program {
	my ($self) = @_;
	
	my $output_format = $self->options('output_format') || 'wig';	
	my $output_file = $self->working_dir() . '/'. $self->job_name.'.'. $output_format;
	
	my @cmd_args;
	
	push @cmd_args, $self->program;
	
	push @cmd_args, '-b', $self->options('background_dir') if ($self->options('background_dir'));
	push @cmd_args, '-c', $self->options('read_count') if ($self->options('read_count'));
	push @cmd_args, '-f', $self->options('fragment_size') if (defined $self->options('fragment_size'));
	push @cmd_args, '-l', $self->options('feature_length') if (defined $self->options('feature_length'));
	push @cmd_args, '-of', $self->options('output_format') if ($self->options('output_format'));
	push @cmd_args, '-p', $self->options('ploidy_dir') if ($self->options('ploidy_dir'));
	push @cmd_args, '-s', $self->options('track_step') if (defined $self->options('track_step'));
	push @cmd_args, '-t', $self->options('threshold') if (defined $self->options('threshold'));
	push @cmd_args, '-v' if ($self->options('verbose'));
	push @cmd_args, '-wg', $self->options('wg_threshold') if ($self->options('wg_threshold'));
	
	my $out_dir = $self->get_temp_dir();
	push @cmd_args, '-o', $out_dir;
		
	push @cmd_args, @{$self->get_input_files_cmd};
	
	my $cmd = join(' ', @cmd_args);
	$self->execute_command_line("bash -c '$cmd'"); # wrapping the cmd like this allows the <(..) file conversion to work
	

	my $merge_cmd = "cat $out_dir/*.$output_format > $output_file";
	
	$self->output_files($output_file);
	$self->bed_file($output_file) if ($output_format eq 'bed');
	
	$self->execute_command_line($merge_cmd);
	
	
}

1;