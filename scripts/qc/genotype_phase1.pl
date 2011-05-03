#!/usr/bin/perl
use strict;
use warnings;


use Data::Dumper;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::GenotypeResultsAdaptor;
use ReseqTrack::Tools::RunAlignment::BFAST;
use ReseqTrack::Tools::RunAlignment::BWA;
use ReseqTrack::Tools::FileSystemUtils qw (create_tmp_process_dir delete_directory);
use ReseqTrack::Tools::QC::GLFTools;
use ReseqTrack::Tools::QC::QCUtils qw (get_params);
use Getopt::Long;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use File::Basename;


my %input;

$input{verbose}        = 0;
$input{working_dir}    = "/nfs/1000g-work/G1K/scratch/rseqpipe/genotype_check/staging/";
$input{subsample_size} = 250000000;
$input{skip_fragment}  = 0;
$input{no_store}       = 0;
my $run_alignment;
my $bams;

$input{update} = 1;

GetOptions (\%input, 'dbhost=s','dbname=s','dbuser=s','dbpass=s',
	    'dbport=s', 'working_dir=s','verbose!','preprocess_exe=s',
	    'type=s','name=s' ,'skip_platform=s', 'debug!',
	    'program=s','snps_bin=s','snps_list=s','reference=s',
	    'validation_method=s','cfg_file=s', 'working_dir=s',
	    'no_store!', 'no_skip_files!','update!',
	   );



if (defined $input{cfg_file}) {
  get_params ($input{cfg_file}, \%input);

}


die "No type given\n"  if (!$input{type});

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
					   -host   => $input{dbhost},
					   -user   => $input{dbuser},
					   -port   => $input{dbport},
					   -dbname => $input{dbname},
					   -pass   => $input{dbpass},
					  );


#  $input{name} (SRR/ERR) is unique in genotype_results table.
my  $gr = $db->get_GenotypeResultsAdaptor;

my $prev_result = $gr->fetch_by_name ($input{name});


if ($prev_result){
  
  throw ( "Already have results for $input{name}. 'update' option not used")
    unless $input{update};
}
else{
  print "No current results in DB\n";
}




my $rmi_a       = $db->get_RunMetaInfoAdaptor;
my $meta_info   = $rmi_a->fetch_by_run_id ($input{name});

throw ("No run_meta_info data found for " .$input{name} . "\n")  if (!$meta_info);

my $paired_length = $meta_info->{paired_length};
my $platform    = $meta_info->{instrument_platform};
print "Platform=$platform\n";

$_= $input{skip_platform};
my @should_skip = /$platform/ig;

if (@should_skip) {
  print "Not running on " . $input{name}." on $platform. Flagged to skip\n";
  sleep (2);
  exit;
}


my $claimed_sample;
$claimed_sample = $meta_info->{sample_name};
print "Claimed sample = $claimed_sample\n";

my $check_snps = have_snps ($claimed_sample, $input{snps_list});


my $ca         = $db->get_CollectionAdaptor;
my $collection = $ca->fetch_by_name_and_type( $input{name}, $input{type} );

if (!$collection) {
  throw "No collection data found for ". $input{name}. " " . $input{type} ."\n";
}

my $other_id = $collection->dbID;
my $files    = $collection->others;


my $BFAST= "ReseqTrack::Tools::RunAlignment::BFAST";

if ($platform eq "ABI_SOLID") {

  $run_alignment = $BFAST->new(
			       -reference      => $input{reference},
			       -program        => $input{program},
			       -working_dir    => $input{working_dir},
			       -samtools       => $input{samtools},
			       -preprocess_exe => $input{preprocess_exe},
			       -subsample_size => $input{subsample_size},
			       -name           => $input{name},
			       -file_info      =>$files,
			       -run_meta_info  =>$meta_info,
			       -input          => $collection,			  
			      );

skip_lane($run_alignment, 
	  $db, 
	  $other_id,
	  $input{name},
	  $claimed_sample,
	  $prev_result,
	 "NO SNPS") if (! $check_snps);

 

 store_skip_short_SOLID($run_alignment, 
			$db, 
			$other_id,
			$input{name},
			$claimed_sample,
		        $prev_result); 


  my $tmp_dir = create_tmp_process_dir ($input{working_dir});

  $run_alignment->working_dir($tmp_dir);

  $run_alignment->decide_file_skip() unless ($input{no_skip_files}) ;
 
  $db->dbc->disconnect_when_inactive(2);

   $run_alignment->subsample_fastq();
  




  my $bang = 0;

  $bang =   $run_alignment->run();
 

  if ($bang > 1) {
    my $tmp_working_dir =  $run_alignment->working_dir();
    delete_directory( $tmp_working_dir);
    throw ("WENT BANG $bang")
  }

  print "Done in bfast genotype section\n";
  my $x = $run_alignment->output_files;
#  print Dumper ($x);
}


##################################

if ( ($platform eq "LS454") || ($platform eq "ILLUMINA")  ) {

  my $BWA="ReseqTrack::Tools::RunAlignment::BWA";

  $run_alignment = $BWA->new(
			     -reference      => $input{reference},
			     -program        => $input{program},
			     -samtools       => $input{samtools},
			     -working_dir    => $input{working_dir},
			     -name           => $input{name},
			     -subsample_size => $input{subsample_size},
			     -file_info      => $files,
			     -paired_length  => $paired_length,
			     -run_meta_info  => $meta_info,
			     -input          => $collection,
			    );

  skip_lane($run_alignment, 
	    $db, 
	    $other_id,
	    $input{name},
	    $claimed_sample,
	    $prev_result,
	    "NO SNPS") if (! $check_snps);


  my $tmp_dir = create_tmp_process_dir ($input{working_dir});

  $run_alignment->working_dir($tmp_dir);

  $run_alignment->decide_file_skip() unless ($input{no_skip_files}) ;

  $db->dbc->disconnect_when_inactive(2);

  $run_alignment->subsample_fastq();
  
  my $bang = 0;
  $bang =   $run_alignment->run();
  
  if ($bang > 1) {
    my $tmp_working_dir =  $run_alignment->working_dir();
    delete_directory( $tmp_working_dir);
    throw ("WENT BANG $bang")
  }
}
##################################


$bams = $run_alignment->output_files;
throw "\n\nNo bams made\n\n" if (! (scalar @$bams) );


if (! $db) {
  $db = ReseqTrack::DBSQL::DBAdaptor->new(
					  -host   => $input{dbhost},
					  -user   => $input{dbuser},
					  -port   => $input{dbport},
					  -dbname => $input{dbname},
					  -pass   => $input{dbpass},
					 );
}




foreach my $bam_file (@$bams){

 
  print $bam_file,"\n";

  my $GLF ="ReseqTrack::Tools::QC::GLFTools";
  my $glf_check = $GLF->new(
			    -reference      => $run_alignment->reference,
			    -program        => $input{glf_program},
			    -name           => $input{name},
			    -verbose        => $input{verbose},   
			    -snps_bin       => $input{snps_bin},
			    -snps_list      => $input{snps_list},
			    -samtools       => $input{samtools},
			    -claimed_sample => $claimed_sample,
			    -working_dir    => $run_alignment->working_dir, 
			    -bam            => $bam_file,
			   );

  $glf_check->run;

  my $results_summary = $glf_check->results_summary;
  throw "no summary" if ( ! $results_summary);



  my $snps    = basename ($input{snps_bin});
  my $aligner = basename ($run_alignment->program);
  my $ref     = basename ($run_alignment->reference);

  my $mapped  = get_bam_stats ($input{samtools},$bam_file);
  $mapped = int ($mapped);




  my $GTR ="ReseqTrack::GenotypeResults";
  my $genotype_results = $GTR->new(
				   -other_id      => $collection->dbID,
				   -table_name    => "collection",
				   -name          => $input{name},
				   -claimed       => $claimed_sample,
				   -top_hit       => $$results_summary{top_hit},
				   -second_hit    => $$results_summary{second_hit},
				   -ratio_2_to_1  => $$results_summary{ratio_2_to_1} ,
				   -ratio_claimed => $$results_summary{ratio_claimed},
				   -reference     => $ref,
				   -snps_bin      => $snps,
				   -aligner       => $aligner,
				   -version       => $run_alignment->program_version,
				   -validation_method =>  $input{validation_method},
				   -percent_mapped => $mapped,
				   -percent_reads_used => $run_alignment->percent_reads_used,
				   -max_bases     => $input{subsample_size},
				   -verdict       => $$results_summary{verdict}, 
				   -cfg_file      => $input{cfg_file},
				   -skip_others   => "1",
				  );






  my  $gra = $db->get_GenotypeResultsAdaptor;

  print Dumper ($genotype_results) if ($input{debug});

  exit if ($input{no_store});

  if ($prev_result){
    $gra->update($genotype_results) ;
  }
  else{
    $gra->store($genotype_results) ;
  }
  
  if ($prev_result){
    $gra->update($genotype_results) ;
  }

}



delete_directory($run_alignment->working_dir) unless ($input{debug});



sub store_skip_short_SOLID{
#TGEN withdrew all 25 and 35mer runs.
  my ( $run_alignment,$db,$other_id,$name,$claimed_sample,$prev_result) = @_;
  my $all_short = 1;

	
  my $read_lengths = $run_alignment->read_lengths;

  foreach my $key (keys %$read_lengths){
    if ( $$read_lengths{$key} > 35){
      $all_short = 0;
    }
  }


  if ($all_short){
    print "Short read lengths for SOLID run.";
    print " Entering dummy entry in genotype_results table.Skipping\n";
  }
  else{
    return;
  }
 

 my $GTR ="ReseqTrack::GenotypeResults";
  my $genotype_results = $GTR->new(
                   -other_id      => $other_id,
                   -table_name    => "collection",
                   -name          => $name,
                   -claimed       => $claimed_sample,
                   -top_hit       => "N/A",
                   -second_hit    => "N/A",
                   -ratio21       => "N/A",
                   -ratio_claimed => "N/A",
                   -reference     => "N/A",
                   -snps_bin      => "N/A",
                   -aligner       => "bfast",
                   -version       => "N/A",
                   -validation_method =>  "N/A",
                   -percent_mapped => 0,
                   -percent_reads_used => 0,
                   -max_bases     => 0,
                   -verdict       => "SHORT_SOLID", 
                   -cfg_file      => "N/A",
                   -skip_others   => "1",
                  );
	
	
  if (! $db) {
    $db = ReseqTrack::DBSQL::DBAdaptor->new(
                       -host   => $input{dbhost},
                       -user   => $input{dbuser},
                      -port   => $input{dbport},
                      -dbname => $input{dbname},
                      -pass   => $input{dbpass},
                     );

	
	}
	
  my  $gra = $db->get_GenotypeResultsAdaptor;

  if ($prev_result){
    $gra->update($genotype_results) ;
  }
  else{
    $gra->store($genotype_results) ;
  }
  
 
  print "Skipping short SOLID run\n";
  exit;
	
}

sub skip_lane{

  my ( $run_alignment,$db,$other_id,$name,$claimed_sample,$prev_result,$reason) = @_;
  
 print "SKIPPING LANE\n";
# print " $run_alignment,$db,$other_id,$name,$claimed_sample,$prev_result,$reason\n";

  
  my $GTR ="ReseqTrack::GenotypeResults";
  my $genotype_results = $GTR->new(
                   -other_id      => $other_id,
                   -table_name    => "collection",
                   -name          => $name,
                   -claimed       => $claimed_sample,
                   -top_hit       => "N/A",
                   -second_hit    => "N/A",
                   -ratio21       => "N/A",
                   -ratio_claimed => "N/A",
                   -reference     => "N/A",
                   -snps_bin      => "N/A",
                   -aligner       => "N/A",
                   -version       => "N/A",
                   -validation_method =>  "N/A",
                   -percent_mapped => 0,
                   -percent_reads_used => 0,
                   -max_bases     => 0,
                   -verdict       => $reason, 
                   -cfg_file      => "N/A",
                   -skip_others   => '1',
                  );
	
	

  if (! $db) {
    $db = ReseqTrack::DBSQL::DBAdaptor->new(
                       -host   => $input{dbhost},
                       -user   => $input{dbuser},
                      -port   => $input{dbport},
                      -dbname => $input{dbname},
                      -pass   => $input{dbpass},
                     );	
	}
	
  my  $gra = $db->get_GenotypeResultsAdaptor;



  if ($prev_result){
    print "Updating\n";
    $gra->update($genotype_results) ;
  }
  else{
    print "Storing\n";
    $gra->store($genotype_results) ;
  }
  
  print "\nSkipping $name:$reason run\n";
  exit;
	
}



sub have_snps{

  my ($sample, $snps_list) = @_;

  if ($sample =~ /unidentified/) {
    print "sample unidentified";
    sleep (2);
    throw ("sample unidentified");
  }

  throw ("Could not find snp list") unless (-e $snps_list);

  my @aa = `grep $sample $snps_list`;
  chomp $aa[0] if @aa;

  if (@aa) {
    print "Have snps for $sample\n";
    return 1;
  }

  print"\nDo not have snps for $claimed_sample\n";

  return 0;
}






sub get_bam_stats {

  my ($samtools, $bam)       = @_;
   
  my $flagstat =  $samtools . " flagstat ";

  my $mapped = 0;

  my $cmd = $flagstat . " " .$bam;
  my @stats = `$cmd`;

 

  foreach (@stats) {
    print;
    if ( $_ =~ /mapped \(/) {
      my @aa = split /\s+/;
      $mapped = $aa[-1];
      $mapped =~ s/\(|\)|\%//g;
    }
  }
        
  return $mapped;
}


 
