package ReseqTrack::Tools::QC::ChkIndelsBAM;

use strict;
use warnings;
use vars qw(@ISA);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Path;
use File::Basename;
use base qw(ReseqTrack::Tools::RunProgram);

=pod

=head1 NAME

ReseqTrack::Tools::QC::ChkIndelsBAM

=head1 SYNOPSIS

Object for running Heng's chk_indel_rg code, which checks BAM files for excessive amount of short indels, which indicates 
low data quality.

=head1 Example

my $chk_indels = ReseqTrack::Tools::QC::ChkIndelsBAM (
                      -input_files => '/path/to/file',
                      -working_dir	=> '/path/to/output_dir
                      );
$chk_indels->run;
my $output_file = $chk_indels->output_files->[0];

=cut


sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);
	return $self;
}

sub run_program {
    my ($self) = @_;
       
    mkpath $self->working_dir unless (-e $self->working_dir);

	print "input file is " . $self->input_files->[0] . "\n";
	my $outfile = $self->working_dir . "/" . basename($self->input_files->[0]) . ".out";
	my $command = $self->program . " " . $self->input_files->[0] . " > " . $outfile;
	
	ReseqTrack::Tools::RunProgram->execute_command_line($command);
	
	$self->output_files($outfile);
	
    return $self;
}

1;