package ReseqTrack::Tools::QC::VerifyBamID;

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;

sub new {
	my ( $class, @args ) = @_;
	my $self = {};
	bless $self, $class;

	my ( $reference, $program, $bam, $claimed_sample, $verbose, $working_dir,
		$plink_file, $bimp_file, $run_meta_info, $samtools)
	  = rearrange(
		[
			qw(
			  REFERENCE
			  PROGRAM
			  BAM
			  CLAIMED_SAMPLE
			  VERBOSE
			  WORKING_DIR
			  PLINK_FILE
                          BIMP_FILE
                          RUN_META_INFO
                          SAMTOOLS
			  )
		],
		@args
	  );

	$self->init();

	$self->bam($bam) if $bam;
	$self->reference($reference);
	$self->program($program);
	$self->claimed_sample($claimed_sample);
	$self->verbose($verbose);
	$self->working_dir($working_dir);
	$self->plink_file($plink_file);
	$self->bimp_file($bimp_file);
        $self->samtools($samtools);
	return $self;
}

sub init {
  my $self = shift;

  $self->reference ("/nfs/1000g-work/G1K/work/bin/VerifyBamID/verifyBamID-0.0.5/reference/human_g1k_v37.fa");
  $self->program ("/nfs/1000g-work/G1K/work/bin/VerifyBamID/verifyBamID-0.0.5/verifyBamID/verifyBamID");
  $self->plink_file("");
  $self->bimp_file ("/nfs/1000g-work/G1K/work/bin/VerifyBamID/verifyBamID-0.0.5/reference/hapmap3_r3_b37_fwd.consensus.HQ.CEU.bimp");

  return;
}

sub run {
  my $self = shift;

#  my $dir = $self->working_dir;
#  chdir($dir)
#    or throw( "Failed to change to " 
#	      . $dir
#	      . " ReseqTrack::Tools::RunAlignment check_dir" );
  $self->construct_run_cmd();
  

  $self->check_bai_present;

  return;

}

######################

sub construct_run_cmd {
	my $self = shift;
	my $cmd;

	$cmd .= $self->program . " ";

	$cmd .= "--reference " . $self->reference . " ";

	$cmd .= "--in " . $self->bam . " ";

	$cmd .= "--bfile " . $self->plink_file   . " " if ($self->plink_file);
	$cmd .= "--bimpfile " . $self->bimp_file . " " if ($self->bimp_file);

        my $out_prefix = $self->working_dir . "\/$$";
	 $out_prefix =~ s/\/\//\//g;
	
	$cmd .= "--out " . $out_prefix . " --verbose";

	print $cmd, "\n";

	print "Running\n";

	eval{
	  `$cmd`;
	};

	if ($@) {
	  throw( "Failed: $cmd" );      
	}
	
	return;
}

######################

sub input_sanity_check {
	my $self = shift;

	throw("Only configured to process single bam or sam")
	  if ( ( defined( $self->bam ) ) && ( defined( $self->sam ) ) );

	throw "No reference sequence specified" if ( !$self->reference );

	throw "No claimed sample id specified" if ( !$self->claimed_sample );

}
######################
sub check_bai_present{
  my $self  = shift;

  my $bam = $self->bam;
  my $bai = $bam . ".bai";

  my $cmd = $self->samtools . " index " . $bam;


  if (! -e $bai){
    print "No bai file present. Creating\n";
    eval{
      `$cmd`;
    };

    if ($@) {
      throw( "Failed: $cmd" );      
    }

  }

 return;
}

sub bam {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{bam} = $arg;
	}
	return $self->{bam};
}

sub reference {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{'reference'} = $arg;
	}
	return $self->{'reference'};
}

sub program {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{'program'} = $arg;
	}
	return $self->{'program'};
}

sub claimed_sample {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{'claimed_sample'} = $arg;
	}
	return $self->{'claimed_sample'};
}

sub verbose {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{'verbose'} = $arg;
	}
	return $self->{'verbose'};
}

sub samtools {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{'samtools'} = $arg;
	}
	return $self->{'samtools'};
}

sub working_dir {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{'working_dir'} = $arg;
	}
	return $self->{'working_dir'};
}

sub plink_file {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{plink_file} = $arg;
	}
	return $self->{plink_file};
}

sub bimp_file {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{bimp_file} = $arg;
	}
	return $self->{bimp_file};
}

sub out_prefix {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{out_prefix} = $arg;
	}
	return $self->{out_prefix};
}

sub run_meta_info {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{run_meta_info} = $arg;
  }
  return $self->{run_meta_info};
}

1;

