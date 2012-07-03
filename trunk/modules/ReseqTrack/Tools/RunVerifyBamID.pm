package ReseqTrack::Tools::RunVerifyBamID;

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use ReseqTrack::Tools::RunProgram;

use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::RunProgram);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $vcf, $out_prefix, $options,$debug,$run_mode ) = rearrange(
		[
		 qw(
			  VCF		 
			  OUT_PREFIX
			  OPTIONS
                          DEBUG
                          RUN_MODE
			  )
		],
		@args
	);

	$self->vcf($vcf);
	$self->out_prefix($out_prefix);
	$self->options($options);
	$self->debug($debug);
	$self->run_mode($run_mode);
	$self->construct_run_cmd();
	$self->add_outfiles();

	return $self;
}

sub run {

	my $self = shift;

	$self->change_dir();
	$self->execute_command_line( $self->command_line );

	return;

}

######################

sub construct_run_cmd {
	my $self = shift;
	my $cmd;
	my $out_prefix;

	my $files = $self->input_files;
	my $bam   = @{$files}[0];
	throw "No bam file name passed\n" if ( !$bam );

	$cmd .= $self->program . " ";

	$cmd .= "--vcf " . $self->vcf . " ";
	$cmd .= "--bam " . $bam . " ";
	$cmd .= $self->options . " " if (defined $self->options );
	$cmd .= "--" . $self->run_mode . " " if ( defined $self->run_mode) ;

	if (defined $self->debug){
	  $cmd .= " --verbose "
	}

	if ( !defined( $self->out_prefix ) ) {
		$out_prefix = $self->working_dir . "\/" . basename($bam);
		$out_prefix =~ s/\.bam//;
		$out_prefix =~ s/\/\//\//g;
	}
	else {
		$out_prefix = $self->out_prefix;
	}

	$cmd .= "--out " . $out_prefix;
	print "\n$cmd\n" 	if (defined $self->debug);
	
	$self->command_line($cmd);
	return;
}

######################
sub command_line {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{command_line} = $arg;
	}
	return $self->{command_line};
}

sub add_outfiles {

	my $self = shift;
	my @output_files;
	my @extensions;

	my $work_dir = $self->working_dir;
	$work_dir =~ s/\/$//;

	my $bam = basename( @{ $self->input_files }[0] );
	$bam =~ s /\.bam$//;
	my $file_base = $work_dir . '/' . $bam;


	@extensions = qw (selfSM selfRG bestSM bestRG );
	
	print "used output files\n";
	foreach my $ext (@extensions) {
		my $filename = $file_base . ".$ext";
		print $self->{$ext} = $filename,"\n";
	}

	return;
}

sub vcf {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{'vcf'} = $arg;
	}
	return $self->{'vcf'};
}

sub out_prefix {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{out_prefix} = $arg;
	}
	return $self->{out_prefix};
}

sub selfSM {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{selfSM} = $arg;
	}
	return $self->{selfSM};
}

sub selfRG {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{selfRG} = $arg;
	}
	return $self->{selfRG};
}

sub bestRG {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{bestRG} = $arg;
	}
	return $self->{bestRG};
}

sub bestSM {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{bestSM} = $arg;
	}
	return $self->{bestSM};
}

sub options {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{options} = $arg;
	}
	return $self->{options};
}

sub debug {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{debug} = $arg;
  }
  return $self->{debug};

}

sub run_mode {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{run_mode} = $arg;
  }
  return $self->{run_mode};

}
1;
