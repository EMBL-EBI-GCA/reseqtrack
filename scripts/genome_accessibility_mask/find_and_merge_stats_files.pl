#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use File::Basename;
use Getopt::Long;
use ReseqTrack::Tools::Loader::File;
use AccessibleGenome::MergeSampleLevelStats;

$| = 1;

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $file_type,
    $base_dir,
    $output_dir,
    $target_count, 
    $chr_list,
    $merge_program,
    $tabix,
    $bgzip,
    $gzip
   );

my $test = 0;
my $run         = 0;
my $verbose     = 0;
my $hostname   = "1000genomes.ebi.ac.uk";

&GetOptions(
  'dbhost=s'    => \$dbhost,
  'dbname=s'    => \$dbname,
  'dbuser=s'    => \$dbuser,
  'dbpass=s'    => \$dbpass,
  'dbport=s'    => \$dbport,
  'base_dir=s'	=> \$base_dir,
  'out_dir=s'	=> \$output_dir,
  'file_type=s'	=> \$file_type,
  'target_count=i'	=> \$target_count,
  'chr_list=s'	=> \$chr_list,
  'run!'                => \$run,
  'verbose!'    		=> \$verbose,
  'hostname=s' 		=> \$hostname,
  'merge_program:s'	=>\$merge_program,
  'tabix:s'		=>\$tabix,
  'bgzip:s'		=>\$bgzip,
  'gzip:s'		=> \$gzip,
  'test!'		=> \$test,
);

$merge_program="/nfs/production/reseq-info/work/bin/genomeMask/statsTools/bin/mergeBaseQCSumStats" unless ($merge_program);
$tabix = "/nfs/production/reseq-info/work/bin/tabix/tabix" unless ($tabix);
$bgzip = "/nfs/production/reseq-info/work/bin/tabix/bgzip" unless ($bgzip);
$gzip = "/usr/bin/gzip" unless ($gzip);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host         => $dbhost,
  -user         => $dbuser,
  -port         => $dbport,
  -dbname       => $dbname,
  -pass         => $dbpass,
);

$db->dbc->disconnect_when_inactive(1);
my $fa = $db->get_FileAdaptor;
my $ca = $db->get_CollectionAdaptor;

my $found_tbi_files = `find $base_dir -name "*.tbi" -print`;

foreach my $found_tbi (split/\n/, $found_tbi_files) {
	#print "$found_tbi\n";
	
	my $stats = $found_tbi;
	$stats =~ s/.tbi//;
	
	if ($fa->fetch_by_name($stats)) {
		#print "file already been loaded $stats\n";
		next;
	}
	else {	
		print "loading $stats .... if run is set\n";
		load_file($stats, $file_type) unless ($test);
	}
}	

my $fos = $fa->fetch_by_type($file_type);
throw("No file object is found; please check file type\n") if (!$fos || @$fos == 0);

my @file_list;
my $test_cnt = 0;
foreach my $fo ( @$fos) {
	#throw ("file " . $fo->name . " not in base dir $base_dir") if ($fo->name !~ /$base_dir/);
	push @file_list, $fo->name;
	$test_cnt++;
	goto TEST if ($test_cnt == 10 && $test);
}

TEST:
my $merger = AccessibleGenome::MergeSampleLevelStats->new(
	-db							=> $db,
	-hostname					=> $hostname,
	-dbhost						=> $dbhost,
	-dbname						=> $dbname,
	-dbuser						=> $dbuser,
	-dbpass						=> $dbpass,
	-dbport						=> $dbport,
	-input_files				=> \@file_list,
	-program					=> $merge_program,
	-tabix						=> $tabix,
	-bgzip						=> $bgzip,
	-gzip						=> $gzip,
	-file_type					=> $file_type,
	-target_count				=> $target_count,
	-chr_list					=> $chr_list,
	-goahead					=> $run,
	-verbose					=> $verbose,
	-working_dir				=> $output_dir,
);

$merger->run;

### SUBS ###
sub load_file {
	my ($file_name, $file_type) = @_;
	
	my $loader = ReseqTrack::Tools::Loader::File->new(						 
						  -file			=> [$file_name],
						  -type			=> $file_type,
						  -hostname		=> $hostname,
						  -dbhost		=> $dbhost,
						  -dbname		=> $dbname,
						  -dbuser		=> $dbuser,
						  -dbpass		=> $dbpass,
						  -dbport		=> $dbport,
						  -do_md5			=> 1,
						  -update_existing	=>0,
						  -die_for_problems	=>1,
						  -assign_type		=>1,
						  #-store_new		=>1, ## this is to force "store" not "update"
						  -debug			=>1,
						  -verbose			=>$verbose,
						 );

	$loader->process_input();
	$loader->create_objects();
	$loader->sanity_check_objects();
	$loader->load_objects() if($run);
	return 1;
}
       
=pod
Synopsis:

The script find_and_merge_stats_files.pl  goes to base_dir, the directory where stats files have been or are being written to by step 1, find all stats files that have been completed, and load them in db with a user defined type (such as "STATS").  Then the script looks for all stats files in the database with the type ("STATS") as input for merge.  It organises them into groups to be merged (group size is defined by --target_count). For each group, it does the following:
                          1. insert a collection in the db
                          2. change the file type to type_ONGOING
                          3. submit a farm job to
                                3a. merge files in the collection
                                3b. change file type to type_MERGED upon successful merge
                                3c. change collection type to type_MERGED upon successful merge
			    3d. load the merged file in file table as original file type, so it can be subject to additional merge

Outputs from find_and_merge_stats_files.pl are merged files recorded in the db with the original file type ("STATS").  The script is run periodically to pick up newly created stats files and newly merged stats files.  Such files will be merged.  The process repeats until one single stats file exist and has nothing to merge to.  At the very end, if the number of files to merge may be smaller than "target_count", you will have to reduce the target count to make the final merge. When the entire db has only one file with the original file type ("STATS"), this one file is the final merge product.

Required parameters:
	-base_dir	the directory to where stats files have been written; this is the same as -final_output_dir in step 1.
	-file_type	newly created stats files by step 1 is stored in db with this file_type ("STATS" is what I used)
	-target_count	the number of files to merge.  As merge stats is memory intensive, I used 5 for the WGS low coverage data.
	-out_dir	where the merged files are written to
	-chr_list	list of chromosomes in the reference genome
	database parameters.

Example:

>perl reseqtrack/scripts/genome_accessibility_mask/find_and_merge_stats_files.pl \
-base_dir /nfs/production/reseq-info/work/zheng/accessible_genome_mask/stats_out_dir_test \
$WRITE_DB_ARGS -dbname zheng_map_1kg_p3_hs38_working \
-file_type TESTSTATS \
-target_count 2 \
-out_dir /nfs/production/reseq-info/work/zheng/accessible_genome_mask/stats_out_dir_test \
-chr_list /nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
-verbose \
-run

Note: If any of the merge job failed (you may check the lsf log files; or when a file with type type_ONGOING stuck for long time, it is a good indication that the file failed), manually remove them from db. 
	a. delete the collection from the collection table;
	b. delete entries from collection_group table
	c. update the file table to change the type from type_ONGOING to type


