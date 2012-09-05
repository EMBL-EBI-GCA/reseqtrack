package ReseqTrack::Tools::QC::FastQScreen;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;

use base qw(ReseqTrack::Tools::RunProgram);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ($keep_text, $keep_graph, $conf_file, $subset, $bowtie_parameters,)	= rearrange( [
	    qw( KEEP_TEXT KEEP_GRAPH CONF_FILE SUBSET BOWTIE_PARAMETERS )
		], @args);

	throw ('FastQScreen configuration file must be specified') unless ($conf_file);
	throw ('FastQScreen configuration file could not be found'.$conf_file) unless (-e $conf_file);

	$self->keep_text($keep_text || 0);
	$self->keep_graph($keep_graph || 0);
	$self->configuration_file($conf_file);
	$self->subset($subset);
	$self->bowtie_parameters($bowtie_parameters);

	return $self;	
}

sub run {
	my ( $self, @args ) = @_;
	
	my $program = $self->program;
	my $working_dir = $self->working_dir;
	my $bowtie_parameters = $self->bowtie_parameters;
	my $subset = $self->subset;
	my $conf_file = $self->configuration_file;
	my @fastq_files = @{$self->input_files};
	
	
	print "Expecting output on ".$self->text_path.$/;
	
	if (-e $self->graph_path){
		unlink $self->graph_path || throw ("Graph output is already present and cannot be deleted - FastQScreen will not run: ".$self->graph_path);
	}
	if (-e $self->text_path){
		unlink $self->text_path || throw ("Text output is already present and cannot be deleted - FastQScreen will not run: ".$self->text_path);
	}
	
	# fastqscreen specifies /usr/bin/perl, which isn't maintained here. Use the currently running perl to execute the script
	my $current_perl_executable = $^X;
	
	my $cmd  = "$current_perl_executable $program --quiet --outdir $working_dir --conf $conf_file ";
	
	$cmd .= '--paired ' if (scalar(@fastq_files) > 1);
	$cmd .= "--subset $subset " if ($subset);
	$cmd .= "--bowtie '$bowtie_parameters' " if ($bowtie_parameters);
		
	$cmd .= join(' ',@fastq_files);
	
	$self->execute_command_line($cmd); 
	
	$self->handle_output(@fastq_files);	
}

sub handle_output {
	my ($self, @fastq_files) = @_;
	
	$self->keep_or_delete($self->keep_graph, $self->graph_path);
	$self->keep_or_delete($self->keep_text, $self->text_path);
}

sub base_output_name {
	my ($self) = @_;
	
	my $working_dir = $self->working_dir;
	my $input_file_name = ${$self->input_files}[0];
	
	my ($name,$path,$suffix) = fileparse($input_file_name);
	
	return $working_dir.'/'.$name.'_screen';
}


sub graph_path {
	my ($self) = @_;
	return $self->base_output_name.'.png';
}

sub text_path {
	my ($self) = @_;
	return $self->base_output_name.'.txt';
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

sub keep_text {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'keep_text'} = $arg;
  }
  return $self->{'keep_text'};
}

sub keep_graph {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'keep_graph'} = $arg;
  }
  return $self->{'keep_graph'};
}

sub configuration_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'configuration_file'} = $arg;
  }
  return $self->{'configuration_file'};
}
sub subset {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'subset'} = $arg;
  }
  return $self->{'subset'};
}
sub bowtie_parameters {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'bowtie_parameters'} = $arg;
  }
  return $self->{'bowtie_parameters'};
}

1;