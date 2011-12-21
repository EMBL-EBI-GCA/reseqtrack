package ReseqTrack::Tools::GATKTools::IndelReAligner;

use strict;
use warnings;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::AlignmentBase;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);

use vars qw(@ISA);
use Data::Dumper;
@ISA = qw(ReseqTrack::Tools::GATKTools);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);
      
	my ($ira_options, $rtc_options,$rtc_knowns,$bam) = rearrange(
		[
		 qw(
                    IRA_OPTIONS
		    RTC_OPTIONS
                    RTC_KNOWNS
                    BAM
			  )
		],
		@args
	);
	$self->bam($bam);

	$self->jar_file ("GenomeAnalysisTK.jar");
     	$self->options ('ira_options' ,$ira_options);
	$self->options ('rtc_options' ,$rtc_options);
	$self->options ('rtc_knowns'  ,$rtc_knowns);

       
	$self->construct_make_intervals_file_cmd();

	if ( ! (defined 	$self->{options}->{'ira_options'})){
	  print "Using default Indel Realigner options\n";
	 	$self->options ('ira_options' ,
				" --LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing ");

	}

	$self->construct_make_indel_realign_cmd();
	

	return $self;
}


sub run {
	my $self = shift;

	$self->create_target_intervals_file();

	$self->create_make_indel_realign_bam();

	return;

}


sub create_make_indel_realign_bam {
  my $self = shift;

  $self->execute_command_line($self->make_indel_realigned_bam_cmd);

  $self->files_to_delete ( $self->bam);

  $self->bam ($self->realigned_bam_name);

  return;
}


sub create_target_intervals_file {
  my $self = shift;

  $self->check_bai_exist($self->bam);

  my $cmd = $self->make_intervals_file_cmd;

  $self->execute_command_line (  $self->make_intervals_file_cmd);
  $self->files_to_delete ( $self->intervals_file);
  
  return;
}


sub construct_make_indel_realign_cmd {
	my $self = shift;

	my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
	$cmd .= $self->gatk_path . "\/" . $self->jar_file;
	$cmd .= " -T IndelRealigner ";

	if ( defined 	$self->{options}->{'ira_options'}) {
	  $cmd .= 	$self->{options}->{'ira_options'};
	}
	else {
	  throw "\nNo RealignerTargetCreator options\n";
	}

	$cmd .= "-targetIntervals " . $self->intervals_file . " ";

	if ( defined 	$self->{options}->{'rtc_knowns'}) {
	  $cmd .= 	$self->{options}->{'rtc_knowns'};
	}
	else {
	  throw "No RealignerTargetCreator knowns\n";
	}

	

	$cmd .= " -R " . $self->reference . " ";
	$cmd .= " -I " . $self->bam;

	my $realigned_bam = $self->bam . ".indel_realigned.bam";
	$cmd .= " -o " . $realigned_bam;
	$self->make_indel_realigned_bam_cmd($cmd);
	print "\n\n$cmd\n\n";

	$self->realigned_bam_name($realigned_bam);

	#$self->files_to_delete($self->bam);	
	#$self->output_bam_files($realigned_bam);

	return;

}
sub construct_make_intervals_file_cmd {
      
	my $self = shift;

	my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
	$cmd .= $self->gatk_path ;
	$cmd .= "\/";
	$cmd .= $self->jar_file;
	$cmd .= " -T RealignerTargetCreator ";

	if ( defined 	$self->{options}->{'rtc_options'}) {
	  $cmd .= 	$self->{options}->{'rtc_options'};
	}
	else {
	  warn "No Make interval files options\n";
	}

	if ( defined 	$self->{options}->{'rtc_knowns'}) {
	  $cmd .= 	$self->{options}->{'rtc_knowns'};
	}
	else {
	  throw "No RealignerTargetCreator knowns\n";
	}



	my $interval_file = $self->bam . ".interval_list ";
	

	$cmd .= " -o $interval_file ";
	$cmd .= " -R " . $self->reference . " ";
	$cmd .= " -I " . $self->bam;

	print "\n$cmd\n";


	$self->intervals_file ($interval_file);
	$self->make_intervals_file_cmd($cmd);
	return;
}

sub intervals_file {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{intervals_file} = $arg;
	}
	return $self->{intervals_file};
}

sub make_intervals_file_cmd {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{make_intervals_file_cmd} = $arg;
	}
	return $self->{make_intervals_file_cmd};
}

sub make_indel_realigned_bam_cmd {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{make_indel_realigned_bam_cmd} = $arg;
	}
	return $self->{make_indel_realigned_bam_cmd};
}


sub realigned_bam_name {
  
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{realigned_bam_name} = $arg;
  }
  return $self->{realigned_bam_name};
}

1;
