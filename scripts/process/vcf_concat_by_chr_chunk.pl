#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use Getopt::Long;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use Time::Local;
use VertRes::Utils::Sam;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::BamUtils;
use File::Path;
use File::Basename;

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $vcf_list,
    $collection_type,
    $collection_name,
    $store,
	$only, #merge vcfs belong to only a certain chr
	$verbose,
);
   
my $output_dir = "./";
my $start_time = time();
my $host_name = "1000genomes.ebi.ac.uk";
my $remote = 0;

&GetOptions( 
  'dbhost=s'      		=> \$dbhost,
  'dbname=s'      		=> \$dbname,
  'dbuser=s'      		=> \$dbuser,
  'dbpass=s'      		=> \$dbpass,
  'dbport=s'      		=> \$dbport,
  'vcf_list:s'			=> \$vcf_list,
  'output_dir:s'		=> \$output_dir,
  'collection_type:s'	=> \$collection_type,
  'collection_name:s'	=> \$collection_name,
  'store!'				=> \$store,
  'only:s'				=> \$only,
  'verbose!'			=> \$verbose,
);    

if ($vcf_list && $collection_type) {
	throw("Please provide either a list of vcfs to concatenate or a collection_type, not both");
}


$output_dir =~ s/\/$//;
mkpath($output_dir) unless (-e $output_dir);
	
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

my $ca = $db->get_CollectionAdaptor;
my $fa = $db->get_FileAdaptor;
my $hist_a = $db->get_HistoryAdaptor;
my $ha = $db->get_HostAdaptor;

my $host = $ha->fetch_by_name($host_name);
if(!$host){
  $host = ReseqTrack::Host->new
      (
       -name => $host_name,
       -remote => $remote
      );
}

$db->dbc->disconnect_when_inactive(1);
	
if ($collection_type && $collection_name)  {
	
	throw("Collection type has to be  %VCF") if ($collection_type !~ /VCF/i);
	throw("Only VCFs based on merged, ALL populations are to be concatenated") if ($collection_name !~ /all/i); # won't concatenate vcfs that are not at all-population level
		
		### FIXME: more work:
		### When the call is made on one super pop a time, a merge is needed between VCFs of different super pop, 
		### each collection is for one chrom chunk and contains VCF files of different super pops
		### files within one collection can be merged to form a file of calls for all pop in this chrom chunk; the merged file should be stored back
		### to collection table as a gatk_all_chr20 type of collection

	if ( $only ) {
		if ( $only =~ /^\d+|^X|^Y|^MT/i) { #if the beginning of the string is digit or X, Y MT
			$only = "chr" . $only;
		}		
		my ($prefix, $algorithm, $sup_pop, $chr) = split (/_/, $collection_name);
		print "chr is $chr, only is $only\n";
		if  ($chr ne $only) {
			throw("Collection $collection_name does not belogn to specified $only");	 
		}
	}	
		
	my $col = $ca->fetch_by_name_and_type($collection_name, $collection_type);
	throw("No collection is found with the name of $collection_name and type $collection_type") if (!$col);
	my $others = $col->others;
	throw("No VCF files found associating with collection " . $col->name) if (!$others || @$others == 0);
	my $vcfs = "";
	my $count = 0;
	foreach my $other ( @$others ) {
		my $vcf_path = $other->name;
		
		throw("VCF file does not exist") unless ( -e $vcf_path );
		$vcfs = $vcfs . " " . $vcf_path;
		print "Input VCF to concat is $vcf_path\n";
		$count++;
	}
	my $out = $output_dir . "/" . $collection_name . ".vcf";
	my $command = "/nfs/1000g-work/G1K/work/bin/vr-codebase/scripts/vcf-concat -s $count $vcfs > $out";
	
	# The -s (sort-merge) option allows the program to open symotaneuosly $count number of files and concat them and sort the resulting file
	# It needs the input vcf files bgzipped and indexed
	 
	print "Running command\n$command\n";
	
	system($command) if ($store);	
	
	my $exit = $?>>8;
	throw("concatinating VCFs failed\n") if ($exit >=1);
	
	`/nfs/1000g-work/G1K/work/bin/tabix/bgzip -f $out`;
	$exit = $?>>8;
	throw("bgzip VCFs failed\n") if ($exit >=1);
	
	my $zipped_file = $out . ".gz";
	my $file_type = $collection_type;
	store_file($zipped_file, $file_type) if ($store);
	
	`/nfs/1000g-work/G1K/work/bin/tabix/tabix -p vcf $zipped_file`;
	$exit = $?>>8;
	throw("bgzip VCFs failed\n") if ($exit >=1);
	
	print "Concatenated file is $zipped_file\n";
		
}		

my $time_elapsed = time() - $start_time;

print "Job took $time_elapsed seconds\n";	

##### SUBS #####
sub store_file {
	my ($file_path, $file_type) = @_;
	my $basename = basename($file_path);
	my $existing_fo = $fa->fetch_by_filename($basename);
	
	if (!$existing_fo || @$existing_fo == 0) {	
		my $fos = create_objects_from_path_list([$file_path], $file_type, $host); #$files is a reference to an array of file objects (in this case, only one element)
		my $fo = $fos->[0]; #to get the first and only file object
	
		if ($store) {
			$fa->store($fo);	
			my $history = ReseqTrack::History->new(
				-other_id 	=> $fo->dbID,
				-table_name => 'file',
				-comment 	=> "concatenated vcfs by chroms or chrom chunks", 
			);
			$hist_a->store($history);
			$fo->history($history);
		}
	}
	else {
		my $new_fo = ReseqTrack::File->new
 	   			(
	      		  -adaptor => $fa,
			      -dbID => $existing_fo->[0]->dbID,
			      -name => $file_path,
#			      -md5 => "NULL",
			      -host => $host,
			      -type => $file_type,
			   );
			
			my $history_ref = $existing_fo->[0]->history;
			my $history;
			
			if (!$history_ref || @$history_ref == 0) {     	  
				$history = ReseqTrack::History->new(
				-other_id => $existing_fo->[0]->dbID,
				-table_name => 'file',
				-comment => "Update concatenated VCF file", 
				);
				$new_fo->history($history);
			}
			else {
				my $comment = calculate_comment($existing_fo->[0], $new_fo); 
				
				if (!$comment) {
					$comment = "fiddling concatenated VCF files, no comment\n";
				}	
				$history = ReseqTrack::History->new(
					-other_id => $existing_fo->[0]->dbID,
					-table_name => 'file',
					-comment => $comment,
				);
				$new_fo->history($history);	
			}
			
			if ($store) {
				$fa->update($new_fo, 1, 1); # the second 1 is allow change name		
			}
	}
	return 1;
}	

=pod

This script merge VCFs by chrom or chrom_chunks. After this merge, there is only one VCF file for each chrom (eg. all_chr1) or one VCF for the whole genome. 
This script runs after running vcf_merge_by_pop to merge vcfs of different populations into one (optional)

This script can be run from collection table, collections with "all" as super population name can be considered; for each collection (all_chr1), the others in the collection_group table are to be merged.

Example:
perl vcf_concat_by_chr_chunk.pl $WRITE_DB_ARGS -dbname zheng_var_call -collection_type DCC_VCF -collection_name samtools_all_chr20 -output_dir /nfs/1000g-work/G1K/work/zheng/snp_calling/samtools -only chr20 -store 


use vcf-concat from vr-codebase
