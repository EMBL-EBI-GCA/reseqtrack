#!/sw/arch/bin/perl 

use strict;

use ReseqTrack::Tools::Bas;  ### Bas module has to be imported first to make sure the enviroment variables are set properly!
use ReseqTrack::Tools::AlignmentBase;

use ReseqTrack::DBSQL::RejectLogAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use File::Basename;
use File::Path;
use Getopt::Long;
use ReseqTrack::Tools::Loader;
use ReseqTrack::Tools::Loader::Archive;
use ReseqTrack::Tools::BamUtils;

#use lib '/nfs/1000g-work/G1K/work/zheng/reseq-personal/zheng/lib';
#use myTIME;
use ReseqTrack::Tools::myTIME;

$| = 1; 

my $start_run = time();

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $help,
	$input_bam_name, #full path of a bam file in archive staging area
   );

my $run 			= 0;
my $verbose 		= 0;
my $out_dir 		= "./";
my $priority 		= 50;
my $host_name 		= "1000genomes.ebi.ac.uk";
my $remote 			= 0;

&GetOptions(
  'out:s'		=> \$out_dir,
  'dbhost=s'    => \$dbhost,
  'dbname=s'    => \$dbname,
  'dbuser=s'    => \$dbuser,
  'dbpass=s'    => \$dbpass,
  'dbport=s'    => \$dbport,
  'help!'		=> \$help,
  'run!' 	 	=> \$run,
  'verbose!'	=> \$verbose,
  'bam=s'		=> \$input_bam_name,
  'priority=i'	=> \$priority,
  'host_name=s'	=> \$host_name,
  'remote!'     => \$remote,
);

if ($help) {
	help_info();
}
if (!$input_bam_name) {
	throw("A bam path in the archive ataging area is required\n");
}

my ($time_stamp, $month_stamp, $day_stamp) = ReseqTrack::Tools::myTIME::get_time();

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host 	=> $dbhost,
  -user 	=> $dbuser,
  -port 	=> $dbport,
  -dbname 	=> $dbname,
  -pass 	=> $dbpass,
);

my $fa = $db->get_FileAdaptor;
my $fo = $fa->fetch_by_name($input_bam_name);

throw("No file object is found; please check the bam file name and type\n") if (!$fo);

my $loga = $db->get_RejectLogAdaptor;

$db->dbc->disconnect_when_inactive(1);

my $bam = $fo->name;
my $file_type = $fo->type;
print "BAM file path is: $bam\n" if ($verbose);

my $bam_basename = basename($bam);
my @name_bits = split (/\./, $bam_basename);
my $ind = $name_bits[0];

###FIXME: fiddle for phase2 BAMs
my $bas_type = "BAS";
if ($file_type eq "NCBI_BAM" )  { ## if it is NCBI BAMs
	$bas_type = "NCBI_BAS";
}
elsif ($file_type eq "EXOME_BAM" ) {
	$bas_type = "EXOME_BAS";
}	
elsif ($file_type eq "EXOME_BI_BAM" ) {
	$bas_type = "EXOME_BI_BAS";
}
elsif (	$file_type eq "EXOME_BC_BAM") {
	$bas_type = "EXOME_BC_BAS";
}
elsif 	($file_type eq "EXOME_BCM_BAM") {
	$bas_type = "EXOME_BCM_BAS";
}
elsif ( $file_type eq "TEST_BAM" ) {
	$bas_type = "TEST_BAS";
}
elsif ($file_type eq "TEST_EXOME_BAM" ) {
	$bas_type = "TEST_EXOME_BAS";
}
elsif ( $file_type eq "P3L_EXOME_BAM" ) {
	$bas_type = "P3L_EXOME_BAS";
}	
elsif ( $file_type eq "P3L_BAM" ) {
	$bas_type = "P3L_BAS";
}	

unless ($bam =~ /\/nfs\/1000g-work\/G1K\/archive_staging\// ) {
	throw("BAM file $bam has to be in archive staging area in order for it to be archived\n");
}

my $ha = $db->get_HostAdaptor;
my $host = $ha->fetch_by_name($host_name);
if(!$host){
  $host = ReseqTrack::Host->new
      (
       -name => $host_name,
       -remote => $remote
      );
}

if ( $file_type !~ /BAM/i ) { ## These is to handle the non-BAM Complete Genomics files
	goto SKIP;
}
	
my $bas_name = $bam . ".bas";
my $bas = $fa->fetch_by_name($bas_name);
my $bas_basename = $bam_basename . ".bas";

######################################################################################################
### Create BAS file if there isn't one in the db and in the dropbox and load it into the database. ###
######################################################################################################
#if (!$bas  ) {  
if (!$bas && $file_type !~ /NCBI/i && $file_type !~ /CG/i ) {  
	my $bas_base = $fa->fetch_by_filename($bas_basename);
	my $found_bas_path;
	
	if ( $bas_base && @$bas_base > 1  ) {
		$found_bas_path = $bas_base->[0]->name;
		write_log($fo, $loga, "FARM: multiple bas files found for $bam. One is $found_bas_path, it is not in the archive_staging area.");
		throw("Multiple bas files found. One is $found_bas_path found for $bam, it is not in archive staging\n");
	}	
	elsif ( $bas_base && @$bas_base == 1 ) {		
		$found_bas_path = $bas_base->[0]->name;
		write_log($fo, $loga, "FARM: an bas file found - $found_bas_path, it is not in the archive_staging area.");
		throw("An bas file found - $found_bas_path, it is not in archive staging area\n");
	}
	else { # create a bas, store in db as new or replace the one exists in db
		my ($sample2, $platform2, $algorithm2, $project2, $analysis2, $chrom2, $date2) = CHECK_AND_PARSE_FILE_NAME($bam);
		#print "chrom is $chrom2\ndate is $date2\nPATH is " . $ENV{"PATH"} . "\nPERL5LIB is " . $ENV{'PERL5LIB'} . "\n";
	                  
		eval{
			my $bas =  ReseqTrack::Tools::Bas->new (
					 -reference =>'/nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/hs37d5/hs37d5.fa',
                     -input_files => $bam,
                     -md5 => $fo->md5,
                     -working_dir => '/nfs/1000g-work/G1K/scratch/zheng/tmp',
                     -need_tags=> 0,
                     -in_parent => 1,
                    );
			#### FIXME: may have to change -need_tags between 0 or 1
			#### older version of reference is -reference =>'/nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/bwa/grc37/human_g1k_v37.fa',
			$bas->run; ## this will create a bas file at the same directory as the bam file, if -in_parent is set to 1; otherwise specify -output_dir
		};
		if($@){
    		print "Cannot generate bas file for $bam\n";
    		print "$@\n";
    		throw("Cannot generate bas file for $bam\n$@");
  		}		
		
		my $hist_a = $db->get_HistoryAdaptor();

		my $md5 = run_md5($bas_name);		
		if (!$bas_base || @$bas_base == 0 ) { ## when there was no old bas object
			my $bas_objs = create_objects_from_path_list([$bas_name], $bas_type, $host); #$files is a reference to an array of file objects (in this case, only one element)
			my $bas_obj = $bas_objs->[0]; #to get the first and only file object
			$bas_obj->md5($md5); #set md5
		 	if ($run) {
				$fa->store($bas_obj);	
				my $history = ReseqTrack::History->new(
					-other_id 	=> $bas_obj->dbID,
					-table_name => 'file',
					-comment 	=> "Created a bas file for BAM that doesn't come with a bas", 
		 		);
		 		$hist_a->store($history);
		 		$bas_obj->history($history);
		 		$bas = $bas_obj;
		 	}
		 	else {
		 	    write_log($fo, $loga, "FARM: -run tag is not set so bas file created is not stored in database.");
				throw("-run tag is not set so bas file created is not stored in database.\n");
		 	}
		}
		else { ## when there is an old bas object
			my $new_bas_object = ReseqTrack::File->new
 	   			(
	      		  -adaptor => $fa,
			      -dbID => $bas_base->[0]->dbID,
			      -name => $bas_name,
			      -md5 => $md5,
			      -host => $host,
			      -type => $bas_type,
			   );
			
			my $history_ref = $bas_base->[0]->history;
			my $history;
			
			if (!$history_ref || @$history_ref == 0) {     	  
				$history = ReseqTrack::History->new(
				-other_id => $bas_base->[0]->dbID,
				-table_name => 'file',
				-comment => "Update bas file", 
				);
				$new_bas_object->history($history);
			}
			else {
				my $comment = calculate_comment($bas_base->[0], $new_bas_object); 
				
				if (!$comment) {
					$comment = "fiddling bas files, no comment\n";
				}	
				$history = ReseqTrack::History->new(
					-other_id => $bas_base->[0]->dbID,
					-table_name => 'file',
					-comment => $comment,
				);
				$new_bas_object->history($history);	
			}
			
			if ($run) {
				$fa->update($new_bas_object, 1, 1); # the second 1 is allow change name		
			}
			else {
				write_log($fo, $loga, "FARM: -run tag is not set so bas file created is not stored in database.");
				throw("-run tag is not set so bas file created is not stored in database.\n");
		 	}		
		}		   
	}	
}

my $bai_name = $bam . ".bai";
my $bai_basename = $bam_basename . ".bai";

my $bai = $fa->fetch_by_name($bai_name);	#this is to query with "name = "$bai_name"
my $bais = $fa->fetch_by_filename($bai_basename); # this is to query with "name like "%$bai_basename"
throw("Have multiple files associated with ".$bai_basename) if (@$bais >= 2);

my $bai_actual = $bais->[0];
my $bai_actual_name = $bai_actual->name;

if(!$bai) {
	
	if ( $bai_actual_name =~ /dropbox/ ) {
		write_log($fo, $loga, "FARM: bai file $bai_actual_name is still in the dropbox, stop processing its bam and bas");
		throw("bai file $bai_actual_name is still in the dropbox, stop processing its bam and bas\n");
	}	
	
	if ($bam =~ /unmapped/ && !$bai_actual) { ## Sanger not uploading .bai files for unmapped bams because they are all size zero with the same md5
		`touch $bai_name`;
		my $hist_a = $db->get_HistoryAdaptor();
		my $bai_objs = create_objects_from_path_list([$bai_name], "BAI", $host); #$files is a reference to an array of file objects (in this case, only one element)
		my $bai_obj = $bai_objs->[0]; #to get the first and only file object
	  	my $md5 = run_md5($bai_name);
		$bai_obj->md5($md5); #set md5
		$fa->store($bai_obj);	
		my $history = ReseqTrack::History->new(
			-other_id 	=> $bai_obj->dbID,
			-table_name => 'file',
			-comment 	=> "Created a bai file for unmapped bams", 
	 	);
	 	$hist_a->store($history);
	 	$bai_obj->history($history);
	 	$bai = $bai_obj;
		#print "bai file generated and stored for $bai_name\n" if $verbose;
	}
	else {
		write_log($fo, $loga, "FARM: No bai file found for mapped BAM file, it is not even in the dropbox");
		throw("No bai file found for mapped BAM file $bam, it is not even in the dropbox\n");
	}	
}	

SKIP:
			
my $archive_list = '/nfs/1000g-work/G1K/scratch/zheng/tmp/' . $bam_basename . ".tmp_archive_list." . $time_stamp;

open (LIST, ">", $archive_list) || throw("Cannot open temparary archive list $archive_list\n");
	
if ( check_this_md5($fo) == 1 ) {
	move_bam_to_trash($db, $fo, $fo->name, $run);
	throw("md5 check failed for $bam, file moved to reject bin\n");
}
elsif ( $bai && check_this_md5($bai) == 1 ) {
	move_bam_to_trash($db, $bai, $bai->name, $run);
	throw("md5 check failed for $bai_name, file moved to reject bin\n");
}	
elsif ( $bas && check_this_md5($bas) == 1  ) { ## No need to do md5check for bas files that have been just created as the md5 was calculated
	move_bam_to_trash($db, $bas, $bas->name, $run);
	throw("md5 check failed for $bas_name, file moved to reject bin\n");
}
else {
	#print "The BAM file $bam and associated bas and bai files passed md5 check\n" if ($verbose);
	if ($bas && $bai) {
		print LIST "$bam\n$bas_name\n$bai_name\n";
	}
	elsif ( $bai ) {
		print LIST "$bam\n$bai_name\n";
	}
	else {
		print LIST "$bam\n";
	}	
	close(LIST);

	my $action_string = "archive";

	my $max_number = 1000;
	
	my $archiver = ReseqTrack::Tools::Loader::Archive->new(
                                                       -list_file => $archive_list,
                                                       -dbhost => $dbhost,
                                                       -dbname => $dbname,
                                                       -dbuser  => $dbuser,
                                                       -dbpass  => $dbpass,
                                                       -dbport  => $dbport,
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
}		

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n" if ($verbose);

`rm $archive_list`; 	

##### SUBS ########
sub help_info {

 exec('perldoc', $0);

}

sub check_this_md5 {
	my ($bam_obj) = @_;
	my $flag = 0;
	my $submitted_md5 = $bam_obj->md5;
	my $bam_file_path = $bam_obj->name;
	
	my $calculated_md5 = run_md5($bam_file_path);  #run_md5 is a FileSystemUtils.pm method
		
	if ($calculated_md5 ne $submitted_md5) {
		write_log($bam_obj, $loga, "FARM: submitted md5: $submitted_md5 and calculated md5: $calculated_md5");
		$flag = 1;
	}
	return $flag;
}

#############################################################################################

=pod

=head1 NAME

 perl -w ~/ReseqTrack/scripts/qc/bam_md5check_and_archive.pl

=head2 Required arguments:

	-dbhost, 			the name of the mysql-host
	-dbname, 			the name of the mysql database
	-dbuser, 			the name of the mysql user
	-dbpass, 			the database password if appropriate
	-dbport, 			the port the mysql instance is running on, this defaults to 4197 the standard
            				 port for mysql-g1kdcc.ebi.ac.uk

	-bam,				a full path of BAM in archive staging area

=head2 Optional arguments:

	-run				when this tag is used, BAM files that have passed md5 check will be archived 
	-out				directory where the log files will be written, default is the run dir
	-priority			archive priority to pass to archive_files.pl; default is 50
	-verbose			default is off, set flag -verbose to print run logs
	-help				this makes the script print out its options
	

=head1 SYNOPSIS

 This script takes one BAM a time, and do the following:
 	- find associated bas and bai file for it
 	- for unmapped BAM, if no bai file exist, create one for it
 	- if no bas file, create bas file
 	- check md5 for each file
 	- archive files that have passed md5 check

=head1 Example:

perl ~/ReseqTrack/scripts/qc/bam_md5check_and_archive.pl -dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -bam /nfs/1000g-work/G1K/archive_staging/test/NA06985/alignment/NA06985.chrom11.ILLUMINA.bwa.SRP0000testArchive.20091216.bam -out somewhere_the_farm_job_writes_log_to -run -priority 99  

perl ~/ReseqTrack/scripts/qc/bam_md5check_and_archive.pl -dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -bam 
/nfs/1000g-work/G1K/archive_staging/test/NA06985/alignment/NA06985.chrom10.ILLUMINA.bwa.CEU.low_coverage.2011b.bam -out /tmp 



