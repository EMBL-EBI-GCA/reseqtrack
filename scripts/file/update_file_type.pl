#!/usr/bin/env perl

#copied fnd edited from ~/reseq-personal/laura/filtered_fastq_0110/scripts/update_type.pl

use strict;
use warnings;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use File::Basename;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;
my $list;
my $input_new_type;
my $new_type_prefix = "WITHDRAWN";
my $run;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'list=s' 		  => \$list,
  'new_type:s' 	  => \$input_new_type, # if new type is WITHDRAWN_BAM,BAS,BAI, no need this tag
  'new_type_prefix:s' 	  => \$new_type_prefix, 
   'run!'		=> \$run,
 );


my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

#throw("Can't update type without a defined type ") unless($new_type);
my $fa = $db->get_FileAdaptor;
my $lines = get_lines_from_file($list);

foreach my $line(@$lines){
  my $name = basename($line);
  #my $files = $fa->fetch_all_like_name($name);
  my $files = $fa->fetch_by_filename($name);
  my $old_file;
  if(!$files || @$files == 0){
    print STDERR $line." doesn't exist in the database\n";
    next;
  }
  if(@$files >= 2){
    print STDERR "Have ".@$files." which match ".$name."\n";
    foreach my $file(@$files){
      print STDERR $file->name."\n";
      if ($file->type =~ /WITHDRAWN/) {
          next;
      }
      else {          
      	$old_file = $file;
      }	
    }
    #next;
  }
  #print $files->[0]->name." ".$files->[0]->type."\n";
  else {
      $old_file = $files->[0];
  }
  
  my $new_file = copy_file_object($old_file);

  my $old_type = $old_file->type;
  
  my $new_type;
 # if ($old_type =~ /PHASE1/) {
 #     $old_type =~ s/PHASE1_//g;
 # }    
      
  if ($input_new_type) {
		$new_file->type($input_new_type);
  }		
  else {	
		$new_type = $new_type_prefix . "_" . $old_type;		
  		$new_file->type($new_type);
  }
  
  $new_file->dbID($old_file->dbID);
  
  print $new_file->name . " new type is $new_type\n"; 
  
  if ( $new_type =~ /WITHDRAWN/ ) {
  	unless($new_file->path =~ /withdrawn/ ){
    	print STDERR "File ".$new_file->path." is not in the withdrawn directory\n";
    	next;
  	}
  }	
  else {
     if($new_file->path =~ /withdrawn/ ){
    	print STDERR "File ".$new_file->path." is in the withdrawn directory while the type is not a WITHDRWAN type\n";
    	next;
  	}
  }	 
  
  my $history = create_history($new_file, $old_file);
  if($history){
    $new_file->history($history);
    $fa->update($new_file) if ($run);
  }else{
    print STDERR "There appears to be no difference between the two file objects\n";
  }
}

=head
perl /homes/zheng/reseq-personal/zheng/bin/update_file_type.pl  -dbhost mysql-g1kdcc -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -list /nfs/1000g-work/G1K/work/zheng/bam_release.20100311.withdraw_list.f2
perl /homes/zheng/reseq-personal/zheng/bin/update_file_type.pl  -dbhost mysql-g1kdcc -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -list /nfs/1000g-work/G1K/work/zheng/bam_release.20100311.withdraw_list.f2.bai -new_type WITHDRAWN_BAI &
