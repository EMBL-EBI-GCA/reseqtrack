#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils;
use Getopt::Long;
use IPC::System::Simple qw(system);
use AccessibleGenome::MergeSampleLevelStats;

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $hostname,
    $merge_program,
    $file_type,
    $file_list,
    $chr_list,    
	$merged_file,
    $col_name,
    $col_type,
	$goahead,
	$verbose,
);


&GetOptions(
  'dbhost=s'    => \$dbhost,
  'dbname=s'    => \$dbname,
  'dbuser=s'    => \$dbuser,
  'dbpass=s'    => \$dbpass,
  'dbport=s'    => \$dbport,
  'hostname=s' 			=> \$hostname,
  'merge_program=s'		=> \$merge_program,
  'file_type=s'			=> \$file_type,
  'file_list=s'			=> \$file_list,
  'merged_file_name=s'	=>\$merged_file,
  'chr_list=s'			=> \$chr_list,
  'col_name=s'			=> \$col_name,
  'col_type=s'			=> \$col_type,
  'goahead!'			=> \$goahead,
  'verbose!'    		=> \$verbose,  
);

my $tabix = "/nfs/production/reseq-info/work/bin/tabix/tabix";
my $bgzip = "/nfs/production/reseq-info/work/bin/tabix/bgzip";
my $gzip = "/usr/bin/gzip";
 
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

my $files = get_lines_from_file($file_list);

my @fos;
foreach my $file (@$files) {
	my $fo = $fa->fetch_by_name($file);
	throw("No file object exists for file $file") if (!$fo);
	throw("File type is " . $fo->type . ", already a merged file") if ($fo->type =~ /MERGED/);
	push @fos, $fo;
}	

my $collection = $ca->fetch_by_name_and_type($col_name, $col_type);

## make a MergeSampleLevelStats obj so I can use the functions in that object
## can pass in minimal number of parameters that are needed for the functions to call 
my $merger = AccessibleGenome::MergeSampleLevelStats->new(
	-dbhost						=> $dbhost,
	-dbname						=> $dbname,
	-dbuser						=> $dbuser,
	-dbpass						=> $dbpass,
	-dbport						=> $dbport,
	-hostname					=> $hostname,	
	-input_files				=> $files,
	-program					=> $merge_program,
	-file_type					=> $file_type,
	-goahead					=> $goahead,
	-verbose					=> $verbose,
);

my $merge_cmd = $merger->program;
$merge_cmd .= " --out $merged_file ";
$merge_cmd .= " --chrList $chr_list ";
$merge_cmd .= join(" ", @$files);

print "merge command is $merge_cmd\n";
system($merge_cmd);

my $exit = $?>>8;  ### these error handlings don't work, as the command itself will throw. 
if ($exit >=1) {
	$merger->update_file_type(\@fos, "", $fa);
	$merger->update_collection_type("ABORTED", $ca);
	throw("merge failed\n"); 
}    

my $bgzip_file = $merged_file . ".bg";
my $bgzip_cmd = "$gzip -cd $merged_file | $bgzip -c >  $bgzip_file";
print "bgzip cmd is $bgzip_cmd\n";
system($bgzip_cmd);
my $exit2 = $?>>8;
if ($exit2 >=1) {
	$merger->update_file_type(\@fos, "", $fa);
	$merger->update_collection_type("ABORTED", $ca);
	throw("bgzip failed\n"); 
}   

my $mv_cmd = "mv $bgzip_file $merged_file";
print "$mv_cmd\n";
system($mv_cmd);
my $exit3 = $?>>8;
if ($exit3 >=1) {
	$merger->update_file_type(\@fos, "", $fa);
	$merger->update_collection_type("ABORTED", $ca);
	throw("mv bgzip failed\n"); 
}   

my $tabix_cmd = "$tabix -s1 -b2 -e2 $merged_file";
print "$tabix_cmd\n";
system($tabix_cmd);
my $exit4 = $?>>8;
if ($exit4 >=1) {
	$merger->update_file_type(\@fos, "", $fa);
	$merger->update_collection_type("ABORTED", $ca);
	throw("tabix failed\n"); 
}   
			
$merger->load_file($merged_file);
$merger->update_file_type(\@fos, "MERGED", $fa);
$merger->update_collection_type("MERGED", $ca, $collection);

unlink($file_list);



  

