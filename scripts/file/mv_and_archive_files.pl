#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use File::Basename;
use File::Path;
use Getopt::Long;
use ReseqTrack::Tools::Loader::File;
use ReseqTrack::Tools::Loader::Archive;
use IPC::System::Simple qw(system);

#### This code only works with the PERL5LIB /nfs/1000g-work/G1K/work/zheng/reseqtrack-code-1059-tags-pre-hive-metadata-merge/modules
#### The Statistics related lines have been commented out in various modules in order for the code to be able to work on the new db schema, 
#### which doesn't have Statistics table

#### This code read a database that supported GRCh38 alignment pipeline, find all crams that have been pushed to sv dropbox,
#### then move them to archive staging area, load in G1K db, check md5 and then archive the ones that passed md5 check.
#### It also changes the file type in the original db to mark the original file can be deleted (_TO_DEL) or should be re-transferred if md5 failed (_TRANSFER_FAILED)

$| = 1;

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $g1k_dbhost,
    $g1k_dbuser,
    $g1k_dbpass,
    $g1k_dbport,
    $g1k_dbname,
    $input_cram_name, #full path of a cram file on EBI server 
   );

my $run         = 0;
my $verbose     = 0;
my $priority    = 50;
my $host_name   = "1000genomes.ebi.ac.uk";
my $move_to_dir = "/nfs/1000g-work/G1K/archive_staging/ftp/data";
my $dropbox_path = "/nfs/1000g-work/G1K/drop/g1k-drop-sv";

&GetOptions(
  'dbhost=s'    => \$dbhost,
  'dbname=s'    => \$dbname,
  'dbuser=s'    => \$dbuser,
  'dbpass=s'    => \$dbpass,
  'dbport=s'    => \$dbport,
  
  'g1k_dbhost=s'    => \$g1k_dbhost,
  'g1k_dbname=s'    => \$g1k_dbname,
  'g1k_dbuser=s'    => \$g1k_dbuser,
  'g1k_dbpass=s'    => \$g1k_dbpass,
  'g1k_dbport=s'    => \$g1k_dbport,
  
  'move_to_dir=s'       => \$move_to_dir,
  'dropbox_path=s'		=> \$dropbox_path,
  'run!'                => \$run,
  'verbose!'    		=> \$verbose,
  'cram=s'              => \$input_cram_name,
  'priority=i'  		=> \$priority,
  'host_name=s' 		=> \$host_name,
);

if (!$input_cram_name) {
        throw("A cram path is required\n");
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host         => $dbhost,
  -user         => $dbuser,
  -port         => $dbport,
  -dbname       => $dbname,
  -pass         => $dbpass,
);

my $g1k_db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host         => $g1k_dbhost,
  -user         => $g1k_dbuser,
  -port         => $g1k_dbport,
  -dbname       => $g1k_dbname,
  -pass         => $g1k_dbpass,
);

$db->dbc->disconnect_when_inactive(1);
$g1k_db->dbc->disconnect_when_inactive(1);

my $fa = $db->get_FileAdaptor;

my $fo = $fa->fetch_by_name($input_cram_name);
throw("No file object is found; please check the cram file name and type\n") if (!$fo);

if ( $fo->type !~ /CRAM/i || $fo->type !~ /PUSHED/i ) { 
        throw("Input file has to be a CRAM file that have been pushed to the dropbox");
}

my $crai_path = $input_cram_name . ".crai";
my $crai_obj = $fa->fetch_by_name($crai_path); 
throw("No file object is found; perhaps the CRAI file hasn't been created yet, $crai_path\n") if (!$crai_obj);

if ( $crai_obj->type !~ /CRAI/i || $crai_obj->type !~ /PUSHED/i ) { 
        throw("This file has to be a CRAI file that have been pushed to the dropbox");
}

move_file_and_load_in_g1k_db($fo);
move_file_and_load_in_g1k_db($crai_obj);

check_md5_change_file_type_and_archive_file($fo);  
check_md5_change_file_type_and_archive_file($crai_obj);

########## SUBS #########
sub move_file_and_load_in_g1k_db {  ## need to delete the original file from ebi nodes, perhaps with a separate script

        my ($file_obj) = @_; #$file is actual file path

        my $filen = basename($file_obj->name);
        my @tmp = split(/\./, $filen);
        my @tmp2 = split(/\_/, $tmp[0]);
        my @tmp3 = split(/-/, $tmp2[0]);
        my $ind = $tmp3[0];
		
	$dropbox_path =~ s/\/$//;
	my $file_path_in_dropbox = $dropbox_path . "/grch38_crams/" . $filen;
        
        $move_to_dir =~ s/\/$//;
        my $new_dir;

        if ($filen =~ /exome/ ) { 
            $new_dir = $move_to_dir . "/" . $ind . "/exome_alignment/";
        }
        elsif ( $filen =~ /low_cov/i )  { 
            $new_dir = $move_to_dir . "/" . $ind . "/alignment/";
        }
        elsif ( $filen =~ /high_cov/i ) {
            $new_dir = $move_to_dir . "/" . $ind . "/high_cov_alignment/";
        }
        else {
            throw("cannot decide directory structure for file $filen");
        }

        unless (-e $new_dir) {
                mkpath($new_dir);
        }

        my $new_file_path = $new_dir . $filen;

	$new_file_path =~ s/\.bam\.cram/\.cram/;

        my $command = "mv $file_path_in_dropbox $new_file_path";    
        print "$command\n" if ($verbose);
        system($command) if ($run);
        my $exit = $?>>8;
        throw("mv failed\n") if ($exit >=1);
		
	my $loader = ReseqTrack::Tools::Loader::File->new(						 
						  -file			=> [$new_file_path],
						  -type			=> $file_obj->type,
						  -hostname		=> $host_name,
						  -dbhost		=> $g1k_dbhost,
						  -dbname		=> $g1k_dbname,
						  -dbuser		=> $g1k_dbuser,
						  -dbpass		=> $g1k_dbpass,
						  -dbport		=> $g1k_dbport,
						  -do_md5			=> 1,
						  -update_existing	=>1,
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

sub check_md5_change_file_type_and_archive_file {
	my ($fo) = @_;
	my $ori_md5 = $fo->md5;
	my $g1k_fos = $g1k_db->get_FileAdaptor->fetch_by_filename(basename($fo->name));
	if ( !$g1k_fos || scalar(@$g1k_fos) == 0) {
		throw("No file is found in the g1k db with the basename of " . basename($fo->name) );
	}	
	my $new_fo;
	foreach my $g1k_fo (@$g1k_fos) {
		print "g1k_type is " . $g1k_fo->type . " ori_type is " . $fo->type . "\n" if ($verbose);
		next if ($g1k_fo->type ne $fo->type);
		$new_fo = $g1k_fo;
	}	
	
	unless ($new_fo->name =~ /\/nfs\/1000g-work\/G1K\/archive_staging\// ) {
        throw("File " . $new_fo->name . " has to be in archive staging area in order for it to be archived\n");
	}
 
	if ($ori_md5 eq $new_fo->md5) {
		my $action_string = "archive";
        my $max_number = 1000;
        my $archiver = ReseqTrack::Tools::Loader::Archive->new(
                                                       -file 	=> [$new_fo->name],
                                                       -dbhost => $g1k_dbhost,
                                                       -dbname => $g1k_dbname,
                                                       -dbuser  => $g1k_dbuser,
                                                       -dbpass  => $g1k_dbpass,
                                                       -dbport  => $g1k_dbport,
                                                       -action => $action_string,
                                                       -verbose => $verbose,
                                                       -priority=>$priority,
                                                       -max_number=>$max_number,
                                                       -no_lock => 1,
                                                      );
        $archiver->process_input();
        #$archiver->cleanup_archive_table($verbose); #leave this out equals to have the -skip_cleanup option
        $archiver->sanity_check_objects();
        $archiver->archive_objects() if $run;
        
        #### Change the file type in the local db so the next process can delete it
        my $new_file_type = $fo->type . "_TO_DEL";
		$fo->type($new_file_type);

        my $history_ref = $fo->history;
        my $history;

        if (!$history_ref || @$history_ref == 0) {
                $history = ReseqTrack::History->new(
                        -other_id => $fo->dbID,
                        -table_name => 'file',
                        -comment => "Move file to a different location",
                );
                $fo->history($history);
        }
        else {
                $history = ReseqTrack::History->new(
                        -other_id => $fo->dbID,
                        -table_name => 'file',
                        -comment => "Change file type",
                );
				$fo->history($history);
        }   
                    		
		$fa->update($fo, 1, 1) if($run); # the second 1 is allow change name
        
	}	
	else {
		
		my $new_file_type = $fo->type . "_TRANSFER_FAILED";
		$fo->type($new_file_type);

        my $history_ref = $fo->history;
        my $history;

        if (!$history_ref || @$history_ref == 0) {
                $history = ReseqTrack::History->new(
                        -other_id => $fo->dbID,
                        -table_name => 'file',
                        -comment => "Move file to a different location",
                );
                $fo->history($history);
        }
        else {
                $history = ReseqTrack::History->new(
                        -other_id => $fo->dbID,
                        -table_name => 'file',
                        -comment => "Change file type",
                );
				$fo->history($history);
        }   
                    		
		$fa->update($fo, 1, 1) if($run); # the second 1 is allow change name
		
		throw("md5 check failed for file " . $new_fo->name);
	}
	return 1;	
}	


=pod

source /nfs/1000g-work/G1K/work/zheng/reseqtrack-code-1059-tags-pre-hive-metadata-merge.sh

perl /nfs/1000g-work/G1K/work/zheng/reseqtrack.20150107/trunk/scripts/file/mv_and_archive_files.pl \
$WRITE_DB_ARGS -dbname zheng_map_1kg_p3_hs38_test \
-g1k_dbhost mysql-g1kdcc-public -g1k_dbuser g1krw -g1k_dbpass thousandgenomes -g1k_dbport 4197 -g1k_dbname g1k_archive_staging_track \
-cram /panfs/nobackup/production/reseq-info/zheng/map_1kg_p3_hs38_by_alt_bwamem_test/IGSR_alignment/NA12347/alignment/NA12347.alt_bwamem_GRCh38DH.20150522.low_cov_test.bam.cram \
-move_to_dir /nfs/1000g-work/G1K/archive_staging/test \
-run \
-verbose


