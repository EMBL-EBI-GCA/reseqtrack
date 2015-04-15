#!/usr/bin/env perl

use strict;
use warnings;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::GeneralUtils;
use File::Path;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Collection;
use Getopt::Long;
use ReseqTrack::Tools::RunSeqtk;
use ReseqTrack::Tools::FileSystemUtils;

$| = 1;

my ($dbhost, $dbuser, $dbpass, $dbname, $dbport);
my $run_id;
my $collection_type;
my $new_collection_type = 'FQ_OK';
my $clobber;
my $output_dir;
my $min_length = 70;
my $program;
my $help;

my $start_run = time();

&GetOptions(
	    'dbhost=s'     => \$dbhost,
	    'dbname=s'     => \$dbname,
	    'dbuser=s'     => \$dbuser,
	    'dbpass=s'     => \$dbpass,
	    'dbport=s'     => \$dbport,
	    'run_id=s' 		=> \$run_id,
	    'output_dir=s' 	=> \$output_dir,
	    'collection_type=s' => \$collection_type,
	    'new_collection_type=s' => \$new_collection_type,
	    'clobber!' 		=> \$clobber,
	    'min_length=s' 	=> \$min_length,
		'program=s'		=> \$program,
	    'help!' 		=> \$help,
	   );

if($help){
  useage();
}

#Creating database connection
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

unless($db){
  throw("Can't run without a database");
}

#Need run id and collection type to get required files and meta info
unless($run_id && $collection_type){
  throw("Can't run without a run id and a collection type");
}

my $ca = $db->get_CollectionAdaptor;
my $collection = $ca->fetch_by_name_and_type($run_id, $collection_type);

throw("Failed to find a collection for ".$run_id." ".$collection_type." from ".$ca->dbc->dbname) 
  unless($collection);
  
#Check if output collection already exists
my $new_collection = $ca->fetch_by_name_and_type($run_id, $new_collection_type);
if($new_collection){
  print $new_collection->name." ".$new_collection_type." already exists\n";
  unless($clobber){
    print "The collection has gone through QA, exiting process\n";
    exit(0);
  }
}

#get input info
my $others = $collection->others;
my $rmia = $db->get_RunMetaInfoAdaptor;
my $rmi = $rmia->fetch_by_run_id($run_id);
throw("Can't run without meta info for ".$run_id) unless($rmi);

#Skip non-ILLUMINA runs
if ($rmi->instrument_platform !~ /ILLUMINA/i) {
	print $run_id." is platform ".$rmi->instrument_platform." so skipping\n";
	$new_collection_type = "NOT_ILLUMINA";
	store_new_collection($new_collection_type);
	exit(0);
}

#Skip non public meta info
if($rmi->status ne 'public'){
  print $run_id." is now ".$rmi->status." in the archive so skipping\n";
  $new_collection_type = "NOT_PUBLIC";
  store_new_collection($new_collection_type);
  exit(0);
}

#Check you have the correct input files
my ($mate1, $mate2, $frag) = assign_files($others);
if($rmi->library_layout eq 'SINGLE'){
  if(!$frag || ($mate1 || $mate2)){
    print "There is a problem for ".$rmi->run_id." is has the wrong files\n";
    print "There are ".@$others." files\n";
    foreach my $file(@$others){
      print $file->name."\n";
    }
    print $rmi->run_id." has library layout ".$rmi->library_layout."\n";
	$new_collection_type = "WRONG_FILE_CNT";
	store_new_collection($new_collection_type);
	exit(0);    
  }
}else{
  unless($mate1 && $mate2){
    print "There is a problem for ".$rmi->run_id." it has the wrong number of files\n";
    foreach my $file(@$others){
      print $file->name."\n";
    }
	$new_collection_type = "WRONG_FILE_CNT";
	store_new_collection($new_collection_type);  
    exit(0);      
  }
####  some paired end ARCHIVE runs each have 3 files - mate1, 2 and frg  SRR393974 
}

if ($mate1 ) {
	throw("File $mate1 doesn't exit") unless (-e $mate1->name);
}

if ( $mate2 ) {
	throw("File $mate2 doesn't exit") unless (-e $mate2->name);
}		

if ( $frag ) {
	throw("File $frag doesn't exit")  unless (-e $frag->name);
}	
 
$db->dbc->disconnect_when_inactive(1);

print "Have ".@$others." from ".$collection->name."\n";

my $read_cnt_for_the_run = $rmi->archive_read_count;
print "From archive: read cnt for the run is $read_cnt_for_the_run\n";

my @inputs;

if($mate1) {
	push @inputs, $mate1->name;
}
if($mate2) {
	push @inputs, $mate2->name;
}
if($frag) {
	push @inputs, $frag->name;
}

#create RunSeqtk object
my $run_seqtk_fqchk = ReseqTrack::Tools::RunSeqtk->new
  (
   -input_files => \@inputs,
   -program => $program,
   -working_dir => $output_dir,
  );
  
$run_seqtk_fqchk->run_fqchk;
my $outs = $run_seqtk_fqchk->output_files;

### seqtk fqchk only count reads that pass the syntax check. So if the read cnt from seqtk is equal to read cnt from archive
### the run passes seqtk fqchk QA.

my $read_cnt_for_the_collection = 0;
foreach my $out ( @$outs ) {
	
	my $lines = get_lines_from_file($out);
	
	my ($min_len_tmp, $max_len_tmp, $avg_len_tmp) = split(/;/, $lines->[0]);
	my ($name1, $min_len) = split(/: /, $min_len_tmp);
	my ($name2, $max_len) = split(/: /, $max_len_tmp);
	my ($name3, $avg_len) = split(/: /, $avg_len_tmp);

	my ($tmp, $total_bases) = split(/\t/, $lines->[2]);
	
	my $total_reads = $total_bases/$avg_len;
	print "By seqtk: read cnt is $total_reads\n";
	$read_cnt_for_the_collection += $total_reads if ($out !~ /\_2\.fastq\.gz/);
	if ($max_len != $min_len) {
		print "Run $run_id failed QA as the reads have different length in $out";
		$new_collection_type = "VAR_READ_LEN";
		store_new_collection($new_collection_type); 
		unlink($out);	 
		exit(0);
	}
	elsif ($avg_len < $min_length) {
		print "Run $run_id failed QA as the reads are shorter than $min_length bp";
		$new_collection_type = "TOO_SHORT";
		store_new_collection($new_collection_type);  
		unlink($out);	
    	exit(0);
	}	
	unlink($out);	
}	
	
print "From seqtk: read cnt for the fastq files in the collection is $read_cnt_for_the_collection\n";

if ($read_cnt_for_the_collection == $read_cnt_for_the_run) {
	print "Run $run_id passed QA\n";
	store_new_collection($new_collection_type);  

}
else {
	print "Run $run_id failed QA as seqtk fqchk read cnt is less than than archive read count";
	$new_collection_type = "FAILED_FQCHK";
	store_new_collection($new_collection_type);  
}	

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";


#### SUBS ####
sub store_new_collection {
	my ($type) = @_;
	my $filtered_collection = ReseqTrack::Collection->new(
							  	  -name => $rmi->run_id,
							      -others => $others,
							      -type => $type,
							      -table_name => 'file',
							     );
	#Store the collection
	$db->get_CollectionAdaptor->store($filtered_collection);
	return 1;
}

		

=pod


perl $ZHENG_RT/scripts/qc/fastq_simple_qa_by_seqtk.pl $WRITE_DB_ARGS -dbname zheng_map_high_cov_reads_to_GRCh38 -run_id SRR826460 -collection_type ARCHIVE_FASTQ -new_collection_type FASTQ_OK -min_length 70 -clobber -output_dir ~/tmp -program /nfs/1000g-work/G1K/work/bin/seqtk/seqtk
 
 
> perl /nfs/production/reseq-info/work/zheng/reseqtrack/trunk/scripts/qc/fastq_simple_qa_by_seqtk.pl $WRITE_DB_ARGS -dbname zheng_map_1kg_p3_hs38 -run_id ERR022062 -collection_type FASTQ -program /nfs/production/reseq-info/work/bin/seqtk/seqtk -output_dir /nfs/gns/homes/zheng/tmp		

-clobber





