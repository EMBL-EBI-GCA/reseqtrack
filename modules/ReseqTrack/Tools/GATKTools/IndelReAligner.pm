=pod

=head1 NAME

ReseqTrack::Tools::GATKTools::IndelReAligner

=head1 SYNOPSIS

Object to create a bam file that in realigned around
known indel sites using GATK RealignerTargetCreator and
IndelRealigner

example

my $REALIGN_AROUND_INDELS = $IR->new(
		     -reference       => $reference,
		     -input_files     => $input{input_files},
		     -rtc_knowns      =>" -known $millsindels.vcf.gz ",
		     -working_dir     => $input{working_dir},
);



=cut

package ReseqTrack::Tools::GATKTools::IndelReAligner;

use strict;
use warnings;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::AlignmentBase;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
#use ReseqTrack::Tools::GATKTools;

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

	my $default_ira =" -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing";

	if ( ! (defined $self->{options}->{'ira_options'})){
	  print "Using default Indel Realigner options\n";
	 	$self->options ('ira_options' ,"$default_ira");
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

  $self->check_bai_exist($self->bam);

  print $self->make_indel_realigned_bam_cmd,"\n";

  $self->execute_command_line($self->make_indel_realigned_bam_cmd);

  $self->files_to_delete ( $self->bam);

  $self->bam ($self->realigned_bam_name);

  return;
}


sub create_target_intervals_file {
  my $self = shift;

  $self->check_bai_exist($self->bam);

  my $cmd = $self->make_intervals_file_cmd;

  print  $self->make_intervals_file_cmd,"\n";

  $self->execute_command_line (  $self->make_intervals_file_cmd);

  $self->files_to_delete ( $self->intervals_file);

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

	#print "\n$cmd\n";

	my $CL = "java jvm_args -jar GenomeAnalysisTK.jar "
	  . "-T RealignerTargetCreator -R \$reference -o \$intervals_file "
	    ."-known \$known_indels_file(s)";

	#For bam header section
	my %PG =('PG'=>'@PG',
		 'ID'=>"gatk_target_interval_creator",
		 'PN'=>"GenomeAnalysisTK",     
		 'PP'=>"sam_to_fixed_bam",
		 'VN'=>"1.2-29-g0acaf2d",
		 'CL'=> $CL,
		 );

	my $COMMENT = '$known_indels_file(s) = '. $self->{options}->{'rtc_knowns'};
	my %CO = ('@CO'=> $COMMENT);

	
	$self->intervals_file ($interval_file);

	print 	$self->intervals_file,"**************\n";
	sleep (5);
	$self->make_intervals_file_cmd($cmd);
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

	$cmd .= " --targetIntervals " . $self->intervals_file . " ";

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
	#print "\n\n$cmd\n\n";

	$self->realigned_bam_name($realigned_bam);

	#$self->files_to_delete($self->bam);
	#$self->output_bam_files($realigned_bam);

	#@PG     ID:bam_realignment_around_known_indels  PN:GenomeAnalysisTK     PP:gatk_target_interval_creator VN:1.2-29-g0acaf2d      CL:java $jvm_args -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference_fasta -I $bam_file -o $realigned_bam_file -targetIntervals $intervals_file -known $known_indels_file(s) -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing




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
