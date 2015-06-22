#!/usr/bin/perl -w

#copied and edited from reseq-personal/laura/filtered_fastq_0110/scripts/update_collection_type.pl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::CollectionUtils;
use File::Basename;
use Getopt::Long;

#use lib '/homes/zheng/reseq-personal/zheng/lib/reseqtrack_hzb/modules/'; #svn check out the ReseqTrack 
use ReseqTrack::Tools::BamUtils;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $list;
my $input_new_type;
my $run = 0;
my $new_type_prefix = "WITHDRAWN";

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'list=s' => \$list,
  'new_type=s' => \$input_new_type,
  'new_type_prefix=s'	=> \$new_type_prefix,
  'run!'		=> \$run,
 );


my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $fa = $db->get_FileAdaptor;

my $lines = get_lines_from_file($list);
my %run_hash; # for fastq
my $to_withdraw_flag = 0;
my %collections_to_withdraw; # for bam

foreach my $line(@$lines){
  chomp $line;
  if ($line =~ /([E|S]RR\d+)/) {
	  $to_withdraw_flag = 1 if ($line =~ /withdrawn/);
	  my $name = basename($line);
	  $name =~ /^([E|S]RR\d+)/;
	  my $run_id = $1;
	  #print $run_id."\n";
	  push(@{$run_hash{$run_id}}, $line);
  }
  elsif ($line =~ /bam/i ) {
      next if ($line =~ /bas/i || $line =~ /bai/i);
      my ($collection_name, $sample, $platform, $algorithm, $project, $analysis, $chrom, $date) = get_collection_name_from_file_name($line);
  	  $collections_to_withdraw{$collection_name} = 1;
  }    
  else {
      warning("What kind of file is this $line\n");
  }    
      	  
}

my $ca = $db->get_CollectionAdaptor;

##### If it is Fastq collections
if (%run_hash) {
	foreach my $run_id(keys(%run_hash)){
	  #print "Have ".$run_id."\n";
	  my $collection;
	  
	  if ($to_withdraw_flag == 1) {
	      $collection = $ca->fetch_by_name_and_type($run_id, 'FILTERED_FASTQ');
	  }
	  else {
	      $collection = $ca->fetch_by_name_and_type($run_id, 'WITHDRAWN_FILTERED_FASTQ');    
	  }
	      
	  if(!$collection){
	    print STDERR "Failed to fetch a collection for ".$run_id.
	        " using FILTERED_FASTQ\n";
	    next;
	  }else{
	    print "Collection " . $collection->name." ".$collection->type."\n";
	    my $total = @{$collection->others};
	    my $diff = 0;
	    my $others = $collection->others;
	    foreach my $other(@{$others}){
	      if($other->type ne $collection->type){
	        $diff++;
	      }
	    }
	    unless($diff == $total){
	      throw("Some of ".$collection->name." types match");
	    }
	    my $new_col = copy_collection_object($collection);
	    if ($input_new_type) {
	        $new_col->type($input_new_type);
	    }
	    else{
		    my $new_type = $others->[0]->type;
		    $new_col->type($new_type);
	    }
	    my $history = create_collection_history_object($new_col, $collection);
	    print $history->comment."\n";
	    if($history){
	      $new_col->history($history) if ($run);
	      $ca->update_type($new_col) if ($run)
	    }
	  }
	}
}
##### If it is BAM collections
elsif ( %collections_to_withdraw ) {
	
	 foreach my $bam_col_name (keys %collections_to_withdraw) { 
	     my $type;
     
	     if($bam_col_name =~ /exome/i ) {
	        if ( $bam_col_name =~ /mosaik/i) {
	         	$type = "EXOME_BC_BAM";
	     	}
	     	else {
	     		$type = "EXOME_BAM";			
	     	}
	     }	  
	     elsif ($bam_col_name =~ /mosaik/i && $bam_col_name !~ /exome/i ) {
	         $type = "NCBI_BAM";
	     }		
	     else {
	         $type = "BAM";
	     }         
	     my $bam_col_obj = $ca->fetch_by_name_and_type($bam_col_name, $type); 
	     if (!$bam_col_obj) {
			warning("No collection object is found for $bam_col_name and $type\n");
	     	next;
	     }
	     print "BAM collection name is $bam_col_name with type $type\n";
		 my $new_col = copy_collection_object($bam_col_obj);
		
		 my $new_type;
		 if ($input_new_type) {
		 	$new_type = $input_new_type;
		 }
		 else {    
		 	$new_type = $new_type_prefix . "_" . $type ;
		 }
		 
		 my $dup_withdrawn_col_obj = $ca->fetch_by_name_and_type($bam_col_name, $new_type);
		 
		 if ( $dup_withdrawn_col_obj ) { # This is to avoid the case when the collection has been withdrawn previously
		     print STDERR "existing collection $bam_col_name with type $new_type\n";
		     $new_type = $new_type . "_2";
		 }    
		 
		 $new_col->type($new_type);		 
		 my $history = create_collection_history_object($new_col, $bam_col_obj);
		 print "Have $bam_col_name ".$history->comment."\n";
		 if($history){
		 	$new_col->history($history) if ($run);
		 	$ca->update_type($new_col) if ($run);
		 }
	 }
}	 	 

=head

perl /homes/zheng/reseq-personal/zheng/bin/update_collection_type.pl -dbhost mysql-g1kdcc -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -list /nfs/1000g-work/G1K/work/zheng/get_fastq_201008/restore_list.archived.f2 -run

EXAMPLE of the input file for the bam:

1000genomes.ebi.ac.uk> head withdrawn_list.files.f1
/nfs/1000g-archive/vol1/ftp/data/NA18690/alignment/NA18690.chrom10.ILLUMINA.bwa.CHD.exon_targetted.20100311.bam
/nfs/1000g-archive/vol1/ftp/data/NA18690/alignment/NA18690.chrom10.ILLUMINA.bwa.CHD.exon_targetted.20100311.bam.bai
/nfs/1000g-archive/vol1/ftp/data/NA18690/alignment/NA18690.chrom10.ILLUMINA.bwa.CHD.exon_targetted.20100311.bam.bas