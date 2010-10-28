package ReseqTrack::Tools::QC::GLFTools;

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

  my (
      $reference,      $program, $snps_list, $snps_bin,
      $bam,            $sam,     $samtools, $glf_file,
      $claimed_sample, $verbose,  $working_dir , $glf_program,
     )
    = rearrange(
		[
		 qw(
     REFERENCE
     PROGRAM
     SNPS_LIST
     SNPS_BIN
     BAM
     SAM
     SAMTOOLS
     GLF_FILE
     CLAIMED_SAMPLE
     VERBOSE
     WORKING_DIR
     GLF_PROGRAM)
		],
		@args
	       );

  $self->bam($bam)           if $bam;
  $self->sam($sam)           if $sam;
  $self->samtools($samtools) if $samtools;
  $self->reference($reference);

  $self->snps_list($snps_list);
  $self->snps_bin($snps_bin);
  $self->program($program);
  $self->glf_program($glf_program);
  $self->glf_file($glf_file);
  $self->claimed_sample($claimed_sample);
  $self->verbose($verbose);
  $self->working_dir($working_dir);


print $working_dir,"\n";
print $self->working_dir,"\n";
  print Dumper ($self);

  return $self;
}
######################
sub run {
  my $self = shift;

  my  $dir = $self->working_dir;
  chdir($dir) or throw("Failed to change to ".$dir.
                       " ReseqTrack::Tools::RunAlignment check_dir");
 

  # just process glfout file
  if ($self->glf_file){
    $self->read_glf_file_to_hash if $self->glf_file;
    $self->evaluate_sample_genotype_results;
    return;
  }


 
  #given bam or sam. Produce/evaluate glf.out results
  if ( $self->sam() || $self->bam() ){
    $self->input_sanity_check();
    $self->process_bam_sam();
    print "glf file ", $self->glf_file, "\n";

    $self->read_glf_file_to_hash if $self->glf_file;

    # process_bam if $self->bam();
    $self->evaluate_sample_genotype_results;
    return;
  }


}
######################
sub run_glf_checkGenotype {
  my $self = shift;

  my $glf     = $self->glf_program();
  my $snps_bin = $self->snps_bin();
  my $glf_file = $self->glf_file;
  my $glfout   = $self->glf_file . ".out";

  my $cmd = 
    $glf . " checkGenotype " . $snps_bin ." " . $glf_file . " > $glfout";
 
  print $cmd,"\n";
#$glf  . " checkGenotype $snps_bin $$.glf  > $glf_out

}
######################
sub check_snps_available{
  my $self = shift;
  my $snps_list  = $self->snps_list;
  my $claimed_sample = $self->claimed_sample;
  
  my @aa = `grep $claimed_sample $snps_list`;
  chomp $aa[0] if (@aa);
  throw "No snps available for $claimed_sample in inventory" if ( !@aa );
 
}
######################

 
######################
sub process_bam_sam {
  my $self = shift;
  print "Process alignment file \n";

  my $file = $self->sam if $self->sam;

  $file = $self->bam if $self->bam;

  my $samtools = $self->samtools;
  my $ref      = $self->reference;

  if ( $self->sam ) {
    print "converting to bam\n";
    my $outbam =  $file;
     $outbam =~ s/sam/bam/;
    my $import_cmd = "$samtools import $ref $file $outbam";
    print $import_cmd, "\n";

    #system()
    $self->bam("$outbam");
  }

  my $bam      = $self->bam();
 

  my $sorted = $self->bam. ".sorted.bam";
  $sorted =~ s/\/\//\//;
# print $bam,"\n";
# print $samtools,"\n";
  my $sort_cmd = "$samtools sort $bam $sorted";
  print $sort_cmd, "\n\n";


  my $index_cmd = "$samtools index $sorted";
  print $index_cmd, "\n";

  
  my $glf_file =  $self->working_dir . '/' . $$ . ".glf";
  my $pileup_cmd = "$samtools pileup -g -f $ref $sorted > $glf_file";
  print $pileup_cmd, "\n";

  
  $self->glf_file("$glf_file");

}
######################
sub evaluate_sample_genotype_results {
  my ($self) = shift;

  my $scores = $self->glf_scores;
  my $sample = $self->claimed_sample;

  my $ratio21       = $$scores{second}{score} / $$scores{top}{score};
  my $claimed_ratio = $$scores{$sample}{score} / $$scores{top}{score};

  if ( $$scores{$sample}{rank} != 1 ) {
#    print "Sample:$sample: is not top ranked:";
#    print "Sample is ranked:",
#      $$scores{$sample}{rank}, "\t", $$scores{$sample}{score};
#    printf "%6.3f\n", $claimed_ratio;
  }

  my $claimed_rank = $$scores{$sample}{rank};

  my $top_sample    = $$scores{top}{sample};
  my $second_sample = $$scores{second}{sample};

 # print "Top Sample $top_sample Top score    :", $$scores{top}{score}, "\n";

  print "RESULTS\:\:\$run_id\:\:";
  print "\:\:1=$top_sample  ";
  print "2=$second_sample ";
  printf "r(2:1)= %3.1f\:\:", $ratio21;
  print "claimed =$sample ";
  printf "r($claimed_rank:1)= %3.1f", $claimed_ratio;

  if ($$scores{$sample}{rank} != 1){
    print "    ::: DOES NOT PASS ";
    print basename ( $self->glf_file() ) ,"\n";
  }
  else{
    print "    :::  PASSES \n";
}
  
  


}
######################
sub read_glf_file_to_hash {
  my ($self) = shift;

  my $file = $self->glf_file();

  my %glf_scores;
  my $rank = 0;

  throw "No glf file" if ( !-e $file );
 # print "Good to go\n";

  open my $IN, '<', $file || die "Nope in glf.out file read";

  while (<$IN>) {

    next if (/entropy/);

    s/,//g;
    s/.snp//g;
    chomp;

    my @a = split /\s+/;

    my $sample = $a[1];
    $sample =~ s/.snp//g;
    my $like  = $a[3];
    my $score = $a[8];

    if ( !defined $glf_scores{$sample} ) {
      $glf_scores{$sample}{like}  = $like;
      $glf_scores{$sample}{score} = $score;
    } else {
      print "Duplicate snp entry for $sample: Ignoring\n";
      next;
    }
    $rank++;

    if ( $rank == 1 ) {
      $glf_scores{top}{score}  = $score;
      $glf_scores{top}{sample} = $sample;
    }
    if ( $rank == 2 ) {
      $glf_scores{second}{score}  = $score;
      $glf_scores{second}{sample} = $sample;
    }

    $glf_scores{$sample}{rank} = $rank;

  }

  if (    $glf_scores{top}{score} <= 0.0001
	  && $glf_scores{second}{score} < 0.0001 ) {
    print "Top score has score of <= 0.0001. Test probably failed\n";
    die "Test failed glf likelihoods <=  0.0001";
  }

  print "Lines read from glf file $rank\n" if $self->verbose;
  $self->glf_scores( \%glf_scores );

}
######################
sub input_sanity_check {
  my $self = shift;

  throw("Only configured to process single bam or sam")
    if ( ( defined( $self->bam ) ) && ( defined( $self->sam ) ) );

  throw "No reference sequence specified" if ( !$self->reference);

  throw "No claimed sample id specified" if ( !$self->claimed_sample );

  throw "No snp list file given" if (! $self->snps_list);

  throw "No snp bin  file given" if (! $self->snps_bin);


}
######################
sub build {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{build} = $arg;
  }
  return $self->{build};
}

sub bam {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{bam} = $arg;
  }
  return $self->{bam};
}

sub sam {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{sam} = $arg;
  }
  return $self->{sam};
}

sub samtools {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{samtools} = $arg;
  }
  return $self->{samtools};
}

sub snps_list {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{snps_list} = $arg;
  }
  return $self->{snps_list};
}

sub snps_bin {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{snps_bin} = $arg;
  }
  return $self->{snps_bin};
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

sub glf_program {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'glf_program'} = $arg;
  }
  return $self->{'glf_program'};
}

sub glf_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'glf_file'} = $arg;
  }
  return $self->{'glf_file'};
}

sub glf_scores {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'glf_scores'} = $arg;
  }
  return $self->{'glf_scores'};
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

sub working_dir {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'working_dir'} = $arg;
  }
  return $self->{'working_dir'};
}


###WRONG ###############


sub set_ncbi36snp_info_from_build {

  my $self = shift;

  my $NCBI36         = "/nfs/1000g-work/G1K/work/bin/ref";
  my $claimed_sample = $self->claimed_sample;
  my $snp_inventory  = "$NCBI36/hapmap_ncbi36_1309.sample_genders";
  my $snp_bin        = "$NCBI36/snps/hapmap3_1309_snps.bin";
  my $ref;

  throw "No:$snp_bin" if (!-e $snp_bin);
  throw "No:$snp_inventory" if (!-e $snp_inventory);



  $self->snp_bin($snp_bin);

  #check if snps exist from this list(1264) . Cannot run without snps.
  $self->snp_inventory($snp_inventory);

  my @a = `grep $claimed_sample $snp_inventory`;
  print "@a\n";
  throw "No snps available for $claimed_sample in NBCI36 inventory" if ( !@a );

  my @b = grep ( /female/i, @a );
  if (@b) {
    $ref = "$NCBI36/human_b36_female.fa";
  } else {
    $ref = "$NCBI36/human_b36_male.fa";
  }

  throw "No:$ref" if (! (-e $ref) );
  $self->reference($ref);

}

#################
sub set_grc37snp_info_from_build {

  my $self = shift;

  my $GRC37         = "/nfs/1000g-work/G1K/work/reference/BWA";
  my $snp_bin       = "$GRC37/GRC37/main_project/snps/snps-v37.june10.bin";
  my $ref           = "$GRC37/GRC37/human_g1k_v37.fa";
  my $snp_inventory = "$GRC37/GRC37/human_g1k_v37.fa";

  throw "No:$snp_bin" if (!-e $snp_bin);
  throw "No:$ref"         if (!-e $ref);
  throw "No:$snp_inventory" if (!-e $snp_inventory);

  $self->snp_bin("$snp_bin");
  $self->reference("$ref");
  $self->snp_inventory("$snp_inventory");

}


1;

