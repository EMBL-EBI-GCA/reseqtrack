#!/usr/bin/perl -w

use strict;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::File;

use File::Basename;
use File::Copy;
use File::Path;

use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;

my $file_list;
my $file_loc;
my $new_root;
my $old_root;
my $output_file;
my $new_type;
my $run = 0;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'file_list=s' => \$file_list,
  'run!' => \$run,
    );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

my $lists = get_lines_from_file($file_list);
  
my $fa = $db->get_FileAdaptor;
	 
PATH:foreach my $list(@$lists){
   my ($path, $new_size) = split /\s+/, $list;
   my $filen = basename($path); 	

   my $file_objs = $fa->fetch_by_filename($filen); 
  
   if(!$file_objs || @$file_objs == 0 ){
     print "Failed to find a file object for ".$filen." in ".$dbname."\n";
     next PATH;
   }
   elsif ( @$file_objs  > 1) {
  	 throw ("two files exist with the same name -  $filen\n");
   }
     
   my $file_obj = $$file_objs[0];
   my $old_size = $file_obj->size;
   my $old_id = $file_obj->dbID;

   if ( $old_size != $new_size ) {
       $file_obj->size($new_size);
       print "changing size from $old_size to $new_size for file " . $file_obj->name . " if (run)\n";
       
       my $history = ReseqTrack::History->new(
			-other_id => $file_obj->dbID,
			-table_name => 'file',
			-comment => "Change file size", 
		 );
		 $file_obj->history($history);
       
       $fa->update($file_obj, 1) if($run);
   }
}   

=head       
  perl ~/reseq-personal/zheng/bin/load_file_size.pl -dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -file_list /nfs/1000g-work/G1K/work/zheng/exome_bam_release_bc_20110521/SOLID/solid.318.samples.bc.md5sum.f2 > foo2 & 
   
