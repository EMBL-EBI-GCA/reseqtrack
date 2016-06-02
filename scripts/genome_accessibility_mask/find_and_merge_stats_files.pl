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

>perl /nfs/production/reseq-info/work/zheng/reseqtrack/scripts/genome_accessibility_mask/find_and_merge_stats_files.pl \
-base_dir /nfs/production/reseq-info/work/zheng/accessible_genome_mask/stats_out_dir_test \
$WRITE_DB_ARGS -dbname zheng_map_1kg_p3_hs38_working \
-file_type TESTSTATS \
-target_count 2 \
-out_dir /nfs/production/reseq-info/work/zheng/accessible_genome_mask/stats_out_dir_test \
-chr_list /nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
-verbose \
