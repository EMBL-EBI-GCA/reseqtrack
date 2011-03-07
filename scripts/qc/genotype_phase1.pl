#!/usr/bin/perl
use strict;
use warnings;


use Data::Dumper;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::GenotypeResultsAdaptor;
use ReseqTrack::Tools::RunAlignment::BFAST;
use ReseqTrack::Tools::RunAlignment::BWA;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::QC::GLFTools;
use ReseqTrack::Tools::QC::QCUtils qw (get_params);
use Getopt::Long;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use File::Basename;


my %input;

$input{verbose}        = 0;
$input{working_dir}    = "/nfs/nobackup/resequencing_informatics/rseqpipe/genotype_check/tmp/staging";
$input{subsample_size} = 250000000;
$input{skip_fragment}  = 0;

my $run_alignment;
my $bams;



GetOptions (\%input, 'dbhost=s','dbname=s','dbuser=s','dbpass=s',
	    'dbport=s', 'working_dir=s','verbose!','preprocess_exe=s',
	    'type=s','name=s' ,'skip_platform=s', 'debug!',
	    'program=s','snps_bin=s','snps_list=s','reference=s',
	    'validation_method=s','cfg_file=s', 'working_dir=s',
	   );



if (defined $input{cfg_file}) {
  get_params ($input{cfg_file}, \%input);

}
#print Dumper (%input);

die "No type given\n"  if (!$input{type});

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
					   -host   => $input{dbhost},
					   -user   => $input{dbuser},
					   -port   => $input{dbport},
					   -dbname => $input{dbname},
					   -pass   => $input{dbpass},
					  );


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

have_snps ($claimed_sample, $input{snps_list});



my $ca         = $db->get_CollectionAdaptor;
my $collection = $ca->fetch_by_name_and_type( $input{name}, $input{type} );

if (!$collection) {
  throw "No collection data found for ". $input{name}. " " . $input{type} ."\n";
}

my $other_id = $collection->dbID;
my $files    = $collection->others;



my  $gr = $db->get_GenotypeResultsAdaptor;
my  $prev_results = $gr->fetch_by_other_id ($other_id);


check_already_run_with_parms ( $prev_results, \%input);


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

  my $tmp_dir = create_tmp_process_dir ($input{working_dir});

  $run_alignment->working_dir($tmp_dir);

  $run_alignment->decide_file_skip();
 
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
  print Dumper ($x);
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

  my $tmp_dir = create_tmp_process_dir ($input{working_dir});

  $run_alignment->working_dir($tmp_dir);

  $run_alignment->decide_file_skip();

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
				  );






  my  $gra = $db->get_GenotypeResultsAdaptor;

  print Dumper ($genotype_results) if ($input{debug});

  $gra->store($genotype_results);
}



delete_directory($run_alignment->working_dir) unless ($input{debug});






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
  sleep (2);
  exit;

  return 0;
}



sub check_already_run_with_parms {

  #do not run if below parameters are found for other_id
  # in genotype_results_table

  my ($geno_results, $curr_parms) = @_;


  my $matches = 0;

  print "Have ", scalar (@$geno_results)," genotype_results\n";

  
  return 0 if ( scalar (@$geno_results) == 0);

  foreach my $result (@$geno_results) {
    $matches = 0;
 
    $matches++ if ($result->max_bases == $$curr_parms{subsample_size} );
    $matches++ if ($result->aligner   eq basename ($$curr_parms{program}) );
    $matches++ if ($result->reference eq basename ($$curr_parms{reference}) );
    $matches++ if (basename ($result->snps_bin)  eq basename ($$curr_parms{snps_bin}) );
    $matches++ if ($result->validation_method    eq $$curr_parms{validation_method} );
    $matches++ if (basename($result->cfg_file) eq basename ($$curr_parms{cfg_file} ));

    print "matches = $matches\n";
    if ($matches == 6) {
      # all params match and test successfully completed 
      print "Found matching gt results\n";
      sleep (2);		# keep farm scheduler happyish
      exit;
    }
  }
  return ;

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


 
