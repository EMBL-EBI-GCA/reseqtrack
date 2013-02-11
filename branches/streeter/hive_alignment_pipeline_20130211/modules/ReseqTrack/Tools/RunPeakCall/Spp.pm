package ReseqTrack::Tools::RunPeakCall::Spp;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunProgram;

use ReseqTrack::Tools::GeneralUtils qw(execute_system_command);

use base qw(ReseqTrack::Tools::RunPeakCall);

=pod

=head1 NAME

ReseqTrack::Tools::RunPeakCall::Spp

=head1 SYNOPSIS

class for running spp. Child class of ReseqTrack::Tools::RunPeakCall

Requires

Takes either a bed file or a bam file.

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
		bin_size => 1,
		separation_range_min => 50,
		separation_range_max => 500,
		window_size => '1e3',
		z_thr => 3,
		verbose => 0,
	};
}

sub new {
	my ( $class, @args ) = @_;
	
	my $self = $class->SUPER::new(@args);

	#setting defaults
	if (!$self->program) {
	   if ($ENV{Rscript}) {
	     $self->program($ENV{R} . '/Rscript');
	   }
	   else {
	     $self->program('Rscript');
	   }
	 }

	return $self;
}

sub can_read_bam {
	return 1;
}

sub can_use_control {
	return 1;
}
sub control_required {
	return 1;
}

sub write_rscript {
	my ($self,$r_script) = @_;
	
	open(my $r_script_fh,'>',$r_script) or throw("Could not open $r_script to write: $!");
	print $r_script_fh $self->r_script_contents();
	close($r_script_fh);
	
	print "Wrote script to $r_script $/ ".$self->r_script_contents().$/ if ($self->options('verbose'));
}
	
sub r_script_contents {
	my ($self) = @_;
	my $bin_size = $self->options('bin_size');
	my $separation_range_min = $self->options('separation_range_min');
	my $separation_range_max = $self->options('separation_range_max');
	my $window_size = $self->options('window_size');
	my $z_thr = $self->options('z_thr');
	
	return <<END;	
args <- commandArgs(trailingOnly = TRUE)
chip.file = args[1]
input.file  = args[2]
output.file	= args[3]

library(spp)
chip.data = read.bam.tags(chip.file)
binding.characteristics <- get.binding.characteristics(chip.data,srange=c($separation_range_min,$separation_range_max),bin=$bin_size);

chip.data <- select.informative.tags(chip.data,binding.characteristics);
chip.data <- remove.local.tag.anomalies(chip.data);

input.data = read.bam.tags(input.file)
input.data <- select.informative.tags(input.data,binding.characteristics);
input.data <- remove.local.tag.anomalies(input.data);

broad.clusters <- get.broad.enrichment.clusters(chip.data,input.data,window.size=$window_size,z.thr=$z_thr,tag.shift=round(binding.characteristics\$peak\$x/2))
write.broadpeak.info(broad.clusters,output.file)
END
}

sub run_program {
	my ($self) = @_;
	
	
	my $temp_dir = $self->get_temp_dir();
	my $r_script = $temp_dir.'/spp.r';
	$self->write_rscript($r_script);	
	
	my $output_file = $self->working_dir() . '/'. $self->job_name.'.spp_bed';
	
	my @cmd_args;
	
	push @cmd_args, $self->program;
	push @cmd_args, $r_script;
	push @cmd_args, @{$self->input_files};
	push @cmd_args, @{$self->control_files};
	push @cmd_args, $output_file;
	
	$self->bed_file($output_file);
	$self->output_files($output_file);
	
	my $cmd = join(' ',@cmd_args);
	#$self->execute_command_line($cmd);
	execute_system_command($cmd);
}

1;