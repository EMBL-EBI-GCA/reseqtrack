#!/usr/bin/env perl
use strict;
use warnings;


use Data::Dumper;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::GenotypeResultsAdaptor;
use ReseqTrack::Tools::RunAlignment::BFAST;
use ReseqTrack::Tools::RunAlignment::BWA;

use ReseqTrack::Tools::FileSystemUtils qw (create_tmp_process_dir delete_directory);

use ReseqTrack::Tools::QC::GLFTools;
use ReseqTrack::Tools::GeneralUtils qw (get_params);
use ReseqTrack::Tools::QC::GLFUtils qw ( check_previous_result  skip_lane
 sample_array_list_of_files have_snps
 decide_file_skip get_max_read_length get_bam_mapped_stats);

use Getopt::Long;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use File::Basename;

my $ALIGNER;
my $run_alignment;
my $bams;
my $prev_result;
my $files_to_delete;

my %input;

$input{verbose}     = 0;
#$input{working_dir} = "/nfs/nobackup/resequencing_informatics/rseqpipe/genotype_check/grc37/staging";
$input{no_store}    = 0;
$input{update} = 1;
$input{save_files_for_deletion} = 0;

GetOptions (\%input, 'dbhost=s','dbname=s','dbuser=s','dbpass=s',
	    'dbport=s', 'working_dir=s','verbose!','preprocess_exe=s',
	    'type=s','name=s' ,'skip_platform=s', 'debug!',
	    'program=s','snps_bin=s','snps_list=s','reference=s',
	    'validation_method=s','cfg_file=s', 'working_dir=s',
	    'no_store!', 'no_skip_files!','update!','subsample_size=s',
	    'file=s@', 'paired_length=s', 'no_subsample!', 'save_files_for_deletion!',
	    'path_must_contain=s',
	   );


if (defined $input{cfg_file}) {
  get_params ($input{cfg_file}, \%input);
}
die "No type given\n"  if (!$input{type});

die "'working_dir' not set\n" if (! (defined $input{working_dir}));


$input{subsample_size} = 250000000 unless (defined $input{subsample_size});



my ($db,$ca,$fa,$gr, $rmi_a ) =  get_db_adaptors (\%input);

$prev_result = check_previous_result ($gr, $input{name}, $input{update}); 

my ($platform,$claimed_sample, $paired_length, $collection) =
  get_run_meta_info ( $rmi_a,$ca, $input{name}, $input{type});


if($platform eq 'COMPLETE_GENOMICS'){
  exit(0);
}

#for event_complete entries. 'skip_platform' listed in cfg file
if (defined $input{skip_platform}){
  $_= $input{skip_platform};
  my @should_skip = /$platform/ig;

  if (@should_skip) {
    print "Not running on " . $input{name}." on $platform. Flagged to skip\n";
    exit;
  }
}



my $have_snps  = have_snps ($claimed_sample, $input{snps_list});
if (!$have_snps) {
  skip_lane ( $db, $collection->dbID , $input{name}, $claimed_sample, $prev_result,
	      "NO SNPS");
  exit;
}


my ($align_these,$need_to_sample) = 
  decide_file_skip ($collection,$input{subsample_size})   unless ($input{no_skip_files}) ;

if (defined $input{path_must_contain}){
  foreach my $check_path(@$align_these){

    if (! ($check_path =~ $input{path_must_contain})){
      my $msg =  "\nfile:$check_path\ndoes not contain required path\n" 
	. $input{path_must_contain} . " <= is required"
	  ."\nCheck all prerequisite events are completed for $input{name}"
	    . " type = " . $input{type};
      
      throw ("$msg\n");
    }
  }
}


my $max_read    = get_max_read_length ($collection, $align_these, 1);


my $tmp_dir     = create_tmp_process_dir ($input{working_dir});


$db->dbc->disconnect_when_inactive(2);


unless ( $input{no_subsample}) {
  if ($need_to_sample){
    ($align_these, $files_to_delete) =
      sample_array_list_of_files( $align_these, $input{subsample_size} , $tmp_dir);
  }
  else{
    print "No need to subsample\n";
  }
}


my %alignment_hash = (
		      -reference      => $input{reference},
		      -program        => $input{program},
		      -samtools       => $input{samtools},
		      -name           => $input{name},
		      -working_dir    => $tmp_dir,
		      -input_files    => $align_these,
		      -save_files_for_deletion =>$input{save_files_for_deletion}, 
		     );

#print Dumper %alignment_hash;




if ($platform eq "ABI_SOLID") {

  $ALIGNER = "ReseqTrack::Tools::RunAlignment::BFAST"; 

 
 
  if ($max_read < 30) {
    skip_lane ( $db, $collection->dbID, $input{name}, $claimed_sample, $prev_result,
		"SHORT_SOLID");

    exit;
  }

  $alignment_hash{-reference}      =  check_reference( $input{reference},$max_read );
  $alignment_hash{-preprocess_exe} =  $input{preprocess_exe};

} elsif ( ($platform eq "LS454") || ($platform eq "ILLUMINA")) {

  $ALIGNER = "ReseqTrack::Tools::RunAlignment::BWA";
  $alignment_hash{-paired_length} =  $paired_length;

}

if ( ! $ALIGNER ) {
  throw  "Do not know what aligner to use for $platform\n";

}


$run_alignment = $ALIGNER->new( %alignment_hash );

$run_alignment->created_files($tmp_dir);

$run_alignment->run();

$bams = $run_alignment->output_files;


throw "\n\nNo bams made\n\n" if (! (scalar @$bams) );


#should only have one bam file.
foreach my $bam_file (@$bams){ 

 

  next if ( ! ( $bam_file =~ /\.bam$/) );

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

  my $mapped  = get_bam_mapped_stats ($input{samtools},$bam_file);
  $mapped     = int ($mapped);


  if (! $db) {
    print "Creating new db adaptor for some reason\n";
    $db = ReseqTrack::DBSQL::DBAdaptor->new(
					    -host   => $input{dbhost},
					    -user   => $input{dbuser},
					    -port   => $input{dbport},
					    -dbname => $input{dbname},
					    -pass   => $input{dbpass},
					   );
  }

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
				   -validation_method =>  $input{validation_method},
				   -percent_mapped => $mapped,
				   -max_bases     => $input{subsample_size},
				   -verdict       => $$results_summary{verdict}, 
				   -cfg_file      => $input{cfg_file},
				   -skip_others   => "1",
				  );

   if ($prev_result ){
    if ($prev_result->verdict ne $genotype_results->verdict){
      dump_comparison ($prev_result, $genotype_results );			
    }		
  }



  my  $gra = $db->get_GenotypeResultsAdaptor;

  print Dumper ($genotype_results) if ($input{debug});

  exit if ($input{no_store});

  if ($prev_result) {
    print "Updating\n";
    $gra->update($genotype_results) ;
  } else {
    print "Storing\n";
    $gra->store($genotype_results) ;
  }
  
  if ($prev_result) {
    $gra->update($genotype_results) ;
  }

}

#delete_directory($run_alignment->working_dir) unless ($input{debug});



##################################


sub get_run_meta_info {
  my ( $rmi_a, $ca, $name,$type ) = @_;
  
  my $meta_info   = $rmi_a->fetch_by_run_id ($name);
  throw ("No run_meta_info data found for " .$name . "\n")  if (!$meta_info);

  my $paired_length = $meta_info->{paired_length};
  my $platform      = $meta_info->{instrument_platform};
  print "$platform  paired length $paired_length\n";
  my $claimed_sample;
  $claimed_sample = $meta_info->{sample_name};
  print "Claimed sample = $claimed_sample\n";

  my $collection = $ca->fetch_by_name_and_type( $name, $type );
  if (!$collection) {
    throw "No collection data found for ". $name. " " . $type ."\n";
  }
  

  return ( $platform,$claimed_sample,$paired_length, $collection);
}



sub check_reference {
  my ($reference, $max_read_length) = @_;
  #  print "Check:",$reference,"\n";
  #  print "Check:",$max_read_length,"\n";
 

  if ( $max_read_length < 40 ) {
    if ( !( $reference =~ /short/i ) ) {
      warning "Reference name looks wrong.No 'short' but read length < 40";
      print "Automatically changing\n";	    
      $reference =~ s/long/short/;
    }
    return $reference
 
  }

  if ( $max_read_length >= 40 ) {
 
    if (!( $reference =~ /long/i)  ) {
      warning "Reference name looks wrong.No 'long' but read length < 40";
      print "Automatically changing\n";
      $reference =~ s/short/long/;
    }
    return $reference;
    
  }


  return  $reference;
}




sub get_db_adaptors {
  my $input  = shift;

  
  my $db = ReseqTrack::DBSQL::DBAdaptor->new(
					     -host   => $$input{dbhost},
					     -user   => $$input{dbuser},
					     -port   => $$input{dbport},
					     -dbname => $$input{dbname},
					     -pass   => $$input{dbpass},
					    );
  my $ca      = $db->get_CollectionAdaptor;
  my $fa      = $db->get_FileAdaptor;
  my $gr      = $db->get_GenotypeResultsAdaptor;
  my $rmi_a   = $db->get_RunMetaInfoAdaptor;

  return ($db, $ca,$fa,$gr,$rmi_a) ;
}


sub dump_comparison{
  my ($prev, $new_results) = @_;

  print "\n================================================\n";
  print "STATUS CHANGE FOR ", $input{name},"\n";
	
  my @cols = qw(top_hit second_hit ratio_2_to_1 ratio_claimed
                percent_mapped verdict);

  print "                  WAS:    NOW\n";      
  foreach my $col ( @cols){
    printf "%15s:",$col;
    print $prev->$col, "\t" ,$new_results->$col,"\n";	            	   
  }
      
  print "\nSNPS used\n";
  print "old:",$prev->snps_bin,"\n";
  print "new:",$new_results->snps_bin,"\n";
	
  print "\n================================================\n";
 return;
}








