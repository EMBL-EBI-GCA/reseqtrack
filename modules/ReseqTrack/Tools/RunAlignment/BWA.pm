package ReseqTrack::Tools::RunAlignment::BWA;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use Data::Dumper;

@ISA = qw(ReseqTrack::Tools::RunAlignment);

=pod

=head1 NAME

=head1 SYNOPSIS

=head1 Example

=head2 new

  Arg [1]   : ReseqTrack:Tools::RunAlignment::BWA
  Arg [2]   : string, options for the single ended bwa mod
  Arg [3]   : string, options for paired end bwa mod
  Function  : create a BWA object, defaults program name to bwa and
  defaults to using ReseqTrack::Tools::RunAlignment::options for sampe options 
  if sampe options aren't defined and options are
  Returntype: ReseqTrack::Tools::RunAlignment::BWA
  Exceptions: n/a
  Example   : my $bwa = ReseqTrack::Tools::RunAlignment::BWA(
                      -program => "bwa",
                      -input => '/path/to/file'
                      -reference => '/path/to/reference',
                      -samtools => '/path/to/samtools',
                      );

=cut

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $samse_options, $sampe_options, $aln_options) =
	  rearrange( [qw(SAMSE_OPTIONS SAMPE_OPTIONS ALN_OPTIONS)], @args );

	#setting defaults
	$self->program('bwa') unless ( $self->program );
	$self->aln_options("-q 15 ") unless ($aln_options);

	
	$self->samse_options($samse_options);
	$self->sampe_options($sampe_options);
	$self->aln_options($aln_options);

	return $self;
}

=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::RunAlignment::BWA
  Arg [2]   : string, options string
  Function  : samse_options and sampe_options are specific bwa commandline
  options and these are the accessor methods for those variables
  Returntype: string
  Exceptions: n/a
  Example   : my $options = $self->samse_options;

=cut


sub run {
    my ($self) = @_;
    $self->change_dir();
    
    if ( $self->fragment_file ) {
        my $sam = $self->run_samse_alignment();
        $self->sam_files($sam);
    }
    
    if ( $self->mate1_file && $self->mate2_file ) {
        my $sam = $self->run_sampe_alignment();
        $self->sam_files($sam);
    }

    $self->run_samtools;

    return;
}

sub run_sampe_alignment {
	my ($self) = @_;
	my $mate1_sai =
	  $self->run_aln_mode( $self->mate1_file, $self->aln_options, "mate1" );
	my $mate2_sai =
	  $self->run_aln_mode( $self->mate2_file, $self->aln_options, "mate2" );
	my $output_sam = $mate1_sai;
	$output_sam =~ s/\.mate1\.sai/\_pe\.sam/;


	if ( $self->paired_length() ) {
		my $paired_length = " -a " . $self->paired_length;
		$self->sampe_options($paired_length);
	}

	my $sampe_cmd =
	    $self->program
	  . " sampe "
	  . $self->sampe_options . " "
	  . $self->reference . " "
	  . $mate1_sai . " "
	  . $mate2_sai . " "
	  . $self->mate1_file . " "
	  . $self->mate2_file . " > "
	  . $output_sam;

	print $sampe_cmd. "\n";
      
        $self->execute_command_line($sampe_cmd);

	$self->files_to_delete($mate1_sai);
	$self->files_to_delete($mate2_sai);

	return $output_sam;
}

sub run_samse_alignment {
	my ($self) = @_;

	my $sai_file =
	  $self->run_aln_mode( $self->fragment_file, $self->aln_options, "frag" );

	my $output_sam = $sai_file;

	$output_sam =~ s/\.frag\.sai/\_se\.sam/;

	my $samse_cmd =
	    $self->program
	  . " samse "
	  . $self->samse_options . " "
	  . $self->reference . " "
	  . $sai_file . " "
	  . $self->fragment_file . " > "
	  . $output_sam;
	  
        $self->execute_command_line($samse_cmd);

	$self->files_to_delete($sai_file);

	return $output_sam;
}

sub run_aln_mode {
	my ( $self, $input_file, $options, $file_ext ) = @_;

	my $bwa_log_file =  $self->working_dir. '/'. $$ .".log";
	$bwa_log_file =~ s/\/\//\//;
	my $do_bwa_log_file = "2>> ". $bwa_log_file;

	$options = $self->aln_options unless ($options);

	my $output_file = $self->working_dir . "/" . $self->job_name;
	$output_file .= "." . $file_ext if ($file_ext);
	$output_file .= ".sai";

	my $aln_command;
	$aln_command .= $self->program;
	$aln_command .= " aln ";
	$aln_command .= $options . "  " if $options;
	$aln_command .= $self->reference . " ";
	$aln_command .= $input_file . " $do_bwa_log_file > ";
	$aln_command .= $output_file;

        $self->execute_command_line($aln_command);

	$self->files_to_delete($bwa_log_file);

	return $output_file;
}

sub samse_options {
    my ( $self, $arg ) = @_;
    unless ( $self->{samse_options} ) {
        $self->{samse_options} = '';
    }
    if ($arg) {
        $self->{samse_options} = $arg;
    }
    return $self->{samse_options};
}

sub sampe_options {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{sampe_options} = $arg;
    }
    return $self->{sampe_options};
}

sub aln_options {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{aln_options} = $arg;
    }
    return $self->{aln_options};
}


1;

