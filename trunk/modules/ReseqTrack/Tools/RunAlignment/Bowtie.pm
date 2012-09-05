package ReseqTrack::Tools::RunAlignment::Bowtie;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use Data::Dumper;
use File::Basename;

@ISA = qw(ReseqTrack::Tools::RunAlignment);

=pod

=head1 NAME

=head1 SYNOPSIS

=head1 Example

=head2 new

  Arg [1]   : ReseqTrack:Tools::RunAlignment::Bowtie
  Arg [2]   : string, options for the single ended bowtie mod
  Arg [3]   : string, options for paired end bowtie mod
  Function  : create a Bowtie object, defaults program name to bowtie and
  defaults to using ReseqTrack::Tools::RunAlignment::options for sampe options 
  if sampe options aren't defined and options are
  Returntype: ReseqTrack::Tools::RunAlignment::BWA
  Exceptions: n/a
  Example   : my $bwa = ReseqTrack::Tools::RunAlignment::BWA(
                      -program => "bowtie",
                      -input => '/path/to/file'
                      -reference => '/path/to/reference',
                      -samtools => '/path/to/samtools',
                      );
=cut

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $options) = rearrange( [qw(OPTIONS)], @args );

	#setting defaults
	$self->program('bowtie') unless ( $self->program );

	$self->options($options);
	
	return $self;
}

=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::RunAlignment::Bowtie
  Arg [2]   : string, options string
  Function  : 
  Returntype: string
  Exceptions: n/a
  Example   : my $options = $self->samse_options;

=cut


sub run {
    my ($self) = @_;
    $self->change_dir();

    if ( $self->fragment_file ) {
		$self->run_samse_alignment;
	}
       
    
    if ( $self->mate1_file && $self->mate2_file ) {
		$self->run_sampe_alignment;
    }

    return;
}

sub run_samse_alignment {
	my ($self) = @_;
	
	my @command = ($self->program);
	push @command, '-S'; #we require SAM format output
	push @command, $self->options if ($self->options);
    push @command, $self->reference;

	if ($self->fragment_file =~ m/\.gz$/){
		unshift @command, 'gzip -dc',$self->fragment_file,'|';
		push @command, '-';
	}
	else {
		push @command, $self->fragment_file;
	}
	
	my $alignment_file = $self->working_dir.'/'.$self->job_name.'.sam';
	
	$self->sam_files($alignment_file);
	
	push @command, $alignment_file;
	
	my $cmd = join ' ',@command;
	$self->execute_command_line($cmd);
	
	
    $self->run_samtools(1);

    return;
}

sub run_sampe_alignment{
	my ($self) = @_;
	
	my @command = ($self->program);
	push @command, '-S'; #we require SAM format output
	push @command, $self->options if ($self->options);
    
	if ($self->paired_length) {
		push @command, '-X';
		push @command, $self->paired_length;
	}
	
	push @command, '--chunkmbs 200';

	push @command, $self->reference;
	
	$self->_add_mate_file(\@command,'-1',$self->mate1_file);
	$self->_add_mate_file(\@command,'-2',$self->mate2_file);

	my $alignment_file = $self->working_dir.'/'.$self->job_name.'.sam';
	
	$self->sam_files($alignment_file);
	
	push @command, $alignment_file;
	
	my $cmd = join ' ',@command;
	#execute in bash to get process substitution
	$self->execute_command_line("bash -c \"$cmd\"");
	
    $self->run_samtools(1);

    return;

}

sub _add_mate_file{
	my ($self,$command_array,$flag,$file) = @_;
	
	if ($file =~ m/\.gz$/){
		my $decompress_cmd = "gzip -dc $file";
		$file = "<($decompress_cmd)";
	}
	
	push @$command_array, $flag, $file;
	
}

sub options {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{options} = $arg;
    }
    return $self->{options};
}
