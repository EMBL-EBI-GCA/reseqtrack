package ReseqTrack::Tools::QC::GLFTools;

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;

sub new {
 my ( $class, @args ) = @_;
 my $self = {};
 bless $self, $class;

 my (
      $reference,      $program, $snp_list, $snp_bin,
      $bam,            $sam,     $samtools, $glf_file,
      $claimed_sample, $verbose, $build,
   )
   = rearrange(
  [
   qw(
     REFERENCE
     PROGRAM
     SNP_LIST
     SNP_BIN
     BAM
     SAM
     SAMTOOLS
     GLF_FILE
     CLAIMED_SAMPLE
     VERBOSE
     BUILD)
  ],
  @args
   );

 $self->bam($bam)           if $bam;
 $self->sam($sam)           if $sam;
 $self->samtools($samtools) if $samtools;
 $self->reference($reference);

 $self->snp_list($snp_list);
 $self->snp_bin($snp_bin);
 $self->program($program);
 $self->glf_file($glf_file);
 $self->claimed_sample($claimed_sample);
 $self->verbose($verbose);
 $self->build($build);

 #print Dumper ($self);
 return $self;
}

sub run {
 my $self = shift;

 $self->input_sanity_check();

 $self->process_bam_sam if $self->sam();

 print "glf file ", $self->glf_file, "\n";

 $self->read_glf_file_to_hash if $self->glf_file;

 # process_bam if $self->bam();
 $self->evaluate_sample_genotype_results;
}

sub run_glf_checkGenotype {
 my $self = shift;

}
sub check_snps_available{
 my $self = shift;
 my $snp_inventory = $self->snp_inventory;
  my $claimed_sample = $self->claimed_sample;
  
  my @a = `grep $claimed_sample $snp_inventory`;

 throw "No snps available for $claimed_sample in inventory" if ( !@a );
 
}

sub set_ncbi36snp_info_from_build {

 my $self = shift;

 my $NCBI36         = "/nfs/1000g-work/G1K/work/bin/ref";
 my $claimed_sample = $self->claimed_sample;
 my $snp_inventory  = "$NCBI36/hapmap_ncbi36_1309.sample_genders";
 my $snp_bin        = "$NCBI36/snps/hapmap3_1309_snps.bin";
 my $ref;

 $self->snp_bin($snp_bin);

 #check if snps exist from this list(1264) . Cannot run without snps.
 $self->snp_inventory($snp_inventory);

 my @a = `grep $claimed_sample $snp_inventory`;

 throw "No snps available for $claimed_sample in NBCI36 inventory" if ( !@a );

 my @b = grep ( /female/i, @a );
 if (@b) {
  $ref = " $NCBI36/human_b36_female.fa";
 }
 else {
  $ref = " $NCBI36/human_b36_male.fa";
 }

 $self->reference($ref);

}
#################
sub set_grc37snp_info_from_build {

 my $self = shift;

 my $GRC37         = "/nfs/1000g-work/G1K/work/reference/BWA";
 my $snp_bin       = "$GRC37/main_project/snps/snps-v37.june10.bin";
 my $ref           = "$GRC37/reference/BWA/GRC37/human_g1k_v37.fa";
 my $snp_inventory = "$GRC37/reference/BWA/GRC37/human_g1k_v37.fa";

throw "No:$snp_bin" if (!-e $snp_bin);
throw "No:$ref"         if (!-e $ref);
throw "No:$snp_inventory" if (!-e $snp_inventory);

 $self->snp_bin($snp_bin);
 $self->reference($ref);
 $self->snp_inventory($snp_inventory);

}

sub process_bam_sam {
 my $self = shift;
 print "Process alignment file \n";

 my $file = $self->sam if $self->sam;
 $file = $self->bam if $self->bam;
 my $samtools = $self->samtools;
 my $ref      = $self->reference;

 if ( $self->sam ) {
  print "converting to bam\n";
  my $import_cmd = "$samtools import $ref $file $$.bam";
  print $import_cmd, "\n";

  #system()
  $self->bam("$$.bam");
 }

 my $bam      = $self->bam();
 my $sort_cmd = "$samtools sort $bam $$.sorted";
 print $sort_cmd, "\n";

 #system ()

 my $index_cmd = "$samtools index $$.sorted.bam";
 print $index_cmd, "\n";

 #system();

 my $pileup_cmd = "$samtools pileup -g -f $ref $$.sorted.bam > $$.glf";
 print $pileup_cmd, "\n";

 #system();
 $self->glf_file("$$.glf");

}

sub evaluate_sample_genotype_results {
 my ($self) = shift;

 my $scores = $self->glf_scores;
 my $sample = $self->claimed_sample;

 my $ratio21       = $$scores{second}{score} / $$scores{top}{score};
 my $claimed_ratio = $$scores{$sample}{score} / $$scores{top}{score};

 if ( $$scores{$sample}{rank} != 1 ) {
  print "Sample:$sample: is not top ranked:";
  print "Sample is ranked:",
    $$scores{$sample}{rank}, "\t", $$scores{$sample}{score};
  printf "%6.3f\n", $claimed_ratio;
 }

 my $claimed_rank = $$scores{$sample}{rank};

 my $top_sample    = $$scores{top}{sample};
 my $second_sample = $$scores{second}{sample};

 print "Top Sample $top_sample Top score    :", $$scores{top}{score}, "\n";

 print "RESULTS\:\:\$run_id\:\:";
 print "\:\:1=$top_sample  ";
 print "2=$second_sample ";
 printf "r(2:1)= %3.1f\:\:", $ratio21;
 print "claimed =$sample ";
 printf "r($claimed_rank:1)= %3.1f\n", $claimed_ratio;

}

sub read_glf_file_to_hash {
 my ($self) = shift;

 my $file = $self->glf_file();

 my %glf_scores;
 my $rank = 0;

 throw "No glf file" if ( !-e $file );
 print "Good to go\n";

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
  }
  else {
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
      && $glf_scores{second}{score} < 0.0001 )
 {
  print "Top score has score of <= 0.0001. Test probably failed\n";
  die "Test failed glf likelihoods <=  0.0001";
 }

 print "Lines read from glf file $rank\n" if $self->verbose;
 $self->glf_scores( \%glf_scores );

}

sub input_sanity_check {
 my $self = shift;

 throw("Only configured to process single bam of sam")
   if ( ( defined( $self->bam ) ) && ( defined( $self->sam ) ) );

 throw "No reference sequence specified" if ( !$self->reference );

 throw "No claimed sample id specified" if ( !$self->claimed_sample );

 throw "Could not figure out which snp.bin file to use"
   if ( !$self->snp_bin && $self->build );

}

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

sub snp_list {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{snp_list} = $arg;
 }
 return $self->{snp_list};
}

sub snp_bin {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{snp_bin} = $arg;
 }
 return $self->{snp_bin};
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

1;

