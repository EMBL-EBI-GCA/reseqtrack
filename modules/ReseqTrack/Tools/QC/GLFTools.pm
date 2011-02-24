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
      $glf_out_file,   $claimed_sample, $verbose,
      $working_dir ,   $name, $skip_sort, $skip_index,
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
     GLF_OUT_FILE
     CLAIMED_SAMPLE
     VERBOSE
     WORKING_DIR
     NAME
     SKIP_SORT
     SKIP_INDEX
     )
		],
		@args
	       );

  $self->name("NONAME");

  $self->bam($bam)           if $bam;
  $self->sam($sam)           if $sam;
  $self->samtools($samtools) if $samtools;
  $self->reference($reference);
  $self->snps_list($snps_list);
  $self->snps_bin($snps_bin);
  $self->program($program);
  $self->glf_out_file($glf_out_file);
  $self->glf_file($glf_file);
  $self->claimed_sample($claimed_sample);
  $self->verbose($verbose);
  $self->working_dir($working_dir);
  $self->name($name);
  $self->skip_sort($skip_sort);
  $self->skip_index($skip_index);

 

  if ( !($self->bam)      &&
       !($self->sam)      &&
       !($self->glf_file) && 
       !($self->glf_out_file) ){

    throw "GLF:Do not know what to do. No input\n";
  }


  return $self;
}


######################
sub run {
  my $self = shift;
  my $have_test_results;
  my $results;
  my $results_summary;

 
  if ($self->glf_out_file) {

    $results = $self->read_to_hash;
    $self->check_glf_run_success ($results);
    $results_summary = $self->sum_results ($results);
    $self->results_summary($results_summary);   
    return 0 ;
  }


  return if (! $self->check_snps_available() );
 

  my  $dir = $self->working_dir;
  chdir($dir) or throw("Failed to change to ".$dir.
                       " ReseqTrack::Tools::RunAlignment check_dir");



  if ($self->glf_file) {

    $self->create_glf_out_file();
    $results = $self->read_to_hash;
    $self->check_glf_run_success ($results);
    $results_summary = $self->sum_results ($results);
    $self->results_summary($results_summary);
    return;
  }

 
  #given bam or sam. Produce/evaluate glf.out results
  if ( $self->sam() || $self->bam() ) {
  
    $self->input_sanity_check();
    $self->process_bam_sam();
    $self->get_bam_stats();
    $self->create_glf_file();
    $self->create_glf_out_file();
 
  }


  if (! $self->glf_out_file) {
    throw ("No glf.out file available ?? Something failed\n");
  }

  $results = $self->read_to_hash;
  $self->check_glf_run_success ($results);
  $results_summary = $self->sum_results ($results);
  $self->results_summary($results_summary);
 
  return 0;
}



#############
sub results_summary {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{results_summary} = $arg;
  }
  return $self->{results_summary};
}


sub sum_results{
  my ($self,$results) = @_;

  my %summary =();
  my $ref = \%summary;

  
  my $top     = parse_glf_lines ($results->{top});
  my $second  = parse_glf_lines ($results->{second});
  my $claimed = parse_glf_lines ($results->{claimed});


  my $r_21       = sprintf "%3.2f", $$second{score}/$$top{score};
  my $r_claimed  = sprintf "%3.2f", $$claimed{score}/$$top{score};

  # print $r_21,"\t", $r_claimed,"\n";
  my $x = $$top{sample};
  $ref->{top_hit}    = $x;
  $ref->{top_sites}    = $$top{sites};
  $ref->{second_hit} = $$second{sample};
  $ref->{ratio21} = $r_21;
  $ref->{ratio_claimed} = $r_claimed;
  $ref->{claimed_sites} = $$claimed{sites};
  $ref->{claimed} = $$claimed{sample};
  
  if ( $$top{sample} eq $$claimed{sample}) {
    $ref->{verdict} = "PASSED";
  } else {
    $ref->{verdict} = "FAILED";
  }

  print "Created glf summary\n";
  return (\%summary);
}
############
sub parse_glf_lines{

  my  ($line) = @_;

  my %glf_data;

  throw "Bad line pass" if (! $line);


  my @aa = (split /\s+/,$line);
  my @bb = split /\./,$aa[1]; #first bit of NA18549.omni.snp or NA18549.snp

  $glf_data{sample} = $bb[0];
  $glf_data{like}   = $aa[3];
  $glf_data{sites}  = $aa[5];
  $glf_data{score}  = $aa[8];


  return (\%glf_data);

}
#############
sub check_glf_run_success {
  my ($self,$results) = @_;

 
  my $top = $$results{top};
  

  my @aa = (split /\s+/,$top);
  my @bb = split /\./,$aa[1]; #first bit of NA18549.omni.snp or NA18549.snp

  my $sample = $bb[0];
  my $like   = $aa[3];
  my $sites  = $aa[5];
  my $score  = $aa[8];



  if ( $score < 0.001 ) {
    die "test failed sites\n";
  }

  if ( $score > 2147483630  ) {
    die "test failed score\n";
  }

 
  return;
}
#############
sub read_to_hash{
  my ($self) = shift;
  my %results;
  my $rank;



  my $file    =  $self->glf_out_file;
  my $claimed =  $self->claimed_sample;



  open my $IN, '<', $file || die "Nope in glf.out file read";

  while (<$IN>) {
  #  print;
    next if (! /sample/);

    s/,//g;

    chomp;

    $rank++;

    print $_,"\n" if ($rank < 9);

    my @aa = split /\s+/;
    my @bb = split /\./,$aa[1];	#first bit of NA18549.omni.snp or NA18549.snp

    my $sample = $bb[0];
 

    $results{top}    = $_ if ($rank   == 1);
    $results{second} = $_ if ($rank   == 2);
    $results{claimed}= $_ if ($sample eq $claimed);
  }


  close ($IN);
  return (\%results);
   
}





############################################
############################################
############################################
sub get_bam_stats {

  my ($self)       = @_;
   
  my $samtools = $self->samtools;
  my $bam      = $self->sorted_bam;


  my $flagstat =  $samtools . " flagstat ";
  my $cmd = $flagstat . " " . $bam;
    
  my @stats = `$cmd`;

  print $bam,"\n";
  foreach (@stats) {
#    print ;
  }
         
  return;
}


sub create_glf_file{

  my ($self) = shift;
  
  my $glf_file =  $self->working_dir . '/' . $self->name . ".glf";
  $glf_file =~ s/\/\//\//;

  my $samtools = $self->samtools();
  my $ref      = $self->reference();
  my $sorted   = $self->sorted_bam();
  my $pileup_cmd = "$samtools pileup -g -f $ref $sorted > $glf_file";
  
  print $pileup_cmd, "\n" if $self->verbose;

  eval {
    `$pileup_cmd`;
  };
  throw "$pileup_cmd failed:$@" if ($@);

  print $glf_file,"\n\n"  if $self->verbose;

  $self->glf_file("$glf_file");

  return;

}

######################
sub create_glf_out_file {
  my $self = shift;

  my $glf      = $self->program();
  my $snps_bin = $self->snps_bin();
  my $glf_file = $self->glf_file;
  my $glfout   = $self->glf_file . ".out";

  my $cmd = 
    $glf . " checkGenotype " . $snps_bin ." " . $glf_file . " > $glfout";
 
  print $cmd,"\n\n";
 
  eval {
    `$cmd`;
  };
 
  throw "$cmd failed:$@" if ($@);
 
  $self->glf_out_file ( $glfout);
  
  return;

}
######################
sub check_snps_available{
  my $self = shift;
  my $snps_list  = $self->snps_list;
  my $claimed_sample = $self->claimed_sample;
  my $name = $self->name;
  #  print "$claimed_sample\n";
  #  print $snps_list,"\n";

  my @aa = `grep $claimed_sample $snps_list`;
  chomp $aa[0] if (@aa);

 
  if ( !@aa ) {

    $self->outcome ("No snps for $claimed_sample $name");

    print "No snps available for $name:$claimed_sample in inventory";
    exit;
  }
  return 1 ;
}


 
######################
sub process_bam_sam {
  my $self = shift;
  my $sorted_bam; 
  my $file;
  print "Process alignment file \n";

  $file = $self->sam if $self->sam;

  $file = $self->bam if $self->bam;

  my $samtools = $self->samtools;
  my $ref      = $self->reference;

  if ( $self->sam ) {
    print $self->sam," ";
    print "Converting to bam\n"  if $self->verbose;

    my $outbam =  $file;
    $outbam =~ s/sam/bam/;
    my $import_cmd = "$samtools import $ref $file $outbam";
    print $import_cmd, "\n" if $self->verbose;

    #system()

    eval {
      `$import_cmd`;
    };
    throw "$import_cmd failed:$@" if ($@);
    
    throw "Failed out make $outbam" if (!-e $outbam);


    $self->bam("$outbam");
  }

  my $bam      = $self->bam();
 
  print "Using BAM $bam\n";

  
  # if you already hav a sorted bam. 
  if ($self->skip_sort) {
    $self->sorted_bam($bam);
    $sorted_bam = $bam;

  } else {
    my $sorted = $self->bam. ".sorted";
    $sorted =~ s/\/\//\//;
    # print $bam,"\n";
    # print $samtools,"\n";
    my $sort_cmd = "$samtools sort $bam $sorted";
    print $sort_cmd, "\n\n"  if $self->verbose;


    eval {
      `$sort_cmd`;
    };
    throw "$sort_cmd failed:$@" if ($@);
    print "============================================\n";
  
    $sorted_bam = $sorted . ".bam";
    $self->sorted_bam($sorted_bam);
  }



  return if ($self->skip_index); # bai file already exists

  my $index_cmd = "$samtools index $sorted_bam";
  print $index_cmd, "\n\n"  if $self->verbose;
 
  eval {
    `$index_cmd`;
  };
  throw "$index_cmd failed:$@" if ($@);

  return;

}


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
  my $res_line;

  $res_line .= "RESULTS\:".  $self->name;
  $res_line .=  "\:1=$top_sample (" . $$scores{top}{sites} . ") ";
  $res_line .=  "2=$second_sample (" . $$scores{second}{sites} . ")";
  my $tmp_r = sprintf " r(2:1)= %3.2f\:", $ratio21;
  $res_line .= $tmp_r;
  $res_line .= " claimed =$sample ";
  $tmp_r = sprintf "r($claimed_rank:1)= %3.2f", $claimed_ratio;
  $res_line .= $tmp_r;

  print $res_line;

  $res_line =~ s/RESULTS://;

  $self->summary($res_line);

  #duplicated runs by Broad have _R suffix.
  if ($top_sample =~ /\_R/) {
    my $rerun =  $top_sample;
    $rerun =~ s/\_R//g;

    if ($rerun eq $sample) {
      $self->outcome("PASSES");
      print " : PASSES Broad R \n";
      return;
    }
  }

  if ($$scores{$sample}{rank} != 1) {
    print " : DOES NOT PASS\n";
    $self->outcome("DOES NOT PASS");
    print basename ( $self->glf_file() ) ,"\n";
  } else {
    $self->outcome("PASSES");
    print " : PASSES \n";
  }
  return ;
}

######################
sub check_against_snp_bin{
  my ($self) = shift;

  my $snp_bin = $self->snps_bin();
  my $glf_out = $self->glf_file . ".out";
  my $run_glf_cmd =  $self->glf_program  . " checkGenotype $snp_bin $$.glf  > $glf_out";
 
  print $run_glf_cmd,"\n" if $self->verbose;

  eval {
    `$run_glf_cmd`;
  };
  throw "$run_glf_cmd failed:$@" if ($@);

  $self->glf_out($glf_out);

  return;
}

######################
sub read_glf_file_to_hash {
  my ($self) = shift;

  my $file = $self->glf_out_file();

  my %glf_scores;
 
  my $rank = 0;


  my $claimed = $self->claimed_sample;
  
  my $name    = $self->name;

  throw "No glf file" if ( !-e $file );


  open my $IN, '<', $file || throw "Nope in glf.out file read";

  while (<$IN>) {

    next if (! /sample/);

    s/,//g;

    chomp;

    my @aa = split /\s+/;
    my @bb = split /\./,$aa[1];	#first bit of NA18549.omni.snp or NA18549.snp

    my $sample = $bb[0];

   
    my $like  = $aa[3];
    my $sites = $aa[5];
    my $score = $aa[8];



    if ( !defined $glf_scores{$sample} ) {
      $glf_scores{$sample}{like}  = $like;
      $glf_scores{$sample}{sites} = $sites;
      $glf_scores{$sample}{score} = $score;
      
     
    } else {
      print "Duplicate snp entry for $sample: Ignoring\n";
      next;
    }

    $rank++;

    if ($rank < 15) {
      print $_,"\n";
    }

    if ($sample eq $claimed) {
      $glf_scores{claimed}{score}  = $score;
      $glf_scores{claimed}{sample} = $sample;
      $glf_scores{claimed}{sites} = $sites;
      $glf_scores{claimed}{line}   = $_;
      $glf_scores{claimed}{rank}  = $rank; 
      
    }


    if ( $rank == 1 ) {
      $glf_scores{top}{sites}  = $sites;
      $glf_scores{top}{score}  = $score;
      $glf_scores{top}{sample} = $sample;
      $glf_scores{top}{line}   = $_;
      
      if ($like == 0.0) {
	print $_;
	$self->summary("FAILED ",$self->name," top likelihood = 0");
	$self->outcome("FAILED");
	print  "\nRESULTS:FAILED Test failed likelihood = 0.0\n";
        return 0;
      }
 
    }
    if ( $rank == 2 ) {
      $glf_scores{second}{sites}  = $sites;  
      $glf_scores{second}{score}  = $score;
      $glf_scores{second}{sample} = $sample;
      $glf_scores{second}{line}   = $_;
    }

    $glf_scores{$sample}{rank} = $rank;

  }
  close ($IN);

  throw ("Failed to read glf.out file") if (! $rank);  

  if ( ($glf_scores{top}{sites} < 1) ||
       ($glf_scores{top}{score} >2147483647) ) {
    $self->summary("FAILED ",$self->name,"  test failed");
    $self->outcome("FAILED");
    print "RESULTS:FAILED ", $self->name," test failed.See $file\n";
    return 0;
  }

  if ( ($glf_scores{top}{sites} == 0) ) {
    $self->summary("FAILED ",$self->name," sites = 0");
    $self->outcome("FAILED");
    print "RESULTS:FAILED", $self->name, " test failed sites = 0";
    return 0;
  }


  if (    $glf_scores{top}{score} <= 0.0001
	  && $glf_scores{second}{score} < 0.0001 ) {
    print "Top score has score of <= 0.0001. Test probably failed\n";

    print  $glf_scores{top}{line};
    print  $glf_scores{second}{line};
    print  "Sample rank = ",  $glf_scores{claimed}{rank},"\n";; 
    print  $glf_scores{claimed}{line};
 
    $self->outcome("FAILED");
    $self->summary("FAILED ",$self->name," glf likelihoods <=  0.0001");
    print "RESULTS:Test failed glf likelihoods <=  0.0001";
    return 0;
  }

  print "Lines read from glf file $rank\n" if $self->verbose;
  $self->glf_scores( \%glf_scores );

  return 1;

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
  return;

}
######################
sub name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{name} = $arg;
  }
  return $self->{name};
}
sub skip_sort {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{skip_sort} = $arg;
  }
  return $self->{skip_sort};
}
sub skip_index {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{skip_index} = $arg;
  }
  return $self->{skip_index};
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

sub sorted_bam {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{sorted_bam} = $arg;
  }
  return $self->{sorted_bam};
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
sub glf_out_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'glf_out_file'} = $arg;
  }
  return $self->{'glf_out_file'};
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


sub result {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'result'} = $arg;
  }
  return $self->{'result'};
}
sub mapped_reads {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'mapped_reads'} = $arg;
  }
  return $self->{'mapped_read'};
}

sub outcome {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'outcome'} = $arg;
  }
  return $self->{'outcome'};
}
sub summary {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'claimed_summary'} = $arg;
  }
  return $self->{'claimed_summary'};
}


1;

