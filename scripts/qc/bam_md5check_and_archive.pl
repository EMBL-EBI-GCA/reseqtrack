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
use ReseqTrack::Tools::GeneralUtils qw(get_time_stamps);

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
	$reference,
	$working_dir,
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
  'reference:s'	=> \$reference,
  'working_dir:s'	=> \$working_dir,
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

my ($time_stamp, $month_stamp, $day_stamp) = get_time_stamps();

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host 	=> $dbhost,
  -user 	=> $dbuser,
  -port 	=> $dbport,
  -dbname 	=> $dbname,
  -pass 	=> $dbpass,
);

my $fa 		= $db->get_FileAdaptor;
my $loga	= $db->get_RejectLogAdaptor;
my $ha		= $db->get_HostAdaptor;

my $host = $ha->fetch_by_name($host_name);
if(!$host){
  $host = ReseqTrack::Host->new
      (
       -name => $host_name,
       -remote => $remote
      );
}

my $fo = $fa->fetch_by_name($input_bam_name);
throw("No file object is found; please check the bam file name and type\n") if (!$fo);

$db->dbc->disconnect_when_inactive(1);

my $bam = $fo->name;
print "BAM file path is: $bam\n" if ($verbose);

unless ($bam =~ /\/nfs\/1000g-work\/G1K\/archive_staging\// ) {
	throw("BAM file $bam has to be in archive staging area in order for it to be archived\n");
}

if ( $fo->type !~ /BAM/i ) { ## These is to handle the non-BAM Complete Genomics files
	goto SKIP;
}

my $bas_name = $bam . ".bas";
my $bas = $fa->fetch_by_name($bas_name);

#######################################################################################################
### Create BAS file if there isn't one in the db and nor in the dropbox; load it into the database. ###
#######################################################################################################
if (!$bas && $fo->type !~ /NCBI/i && $fo->type !~ /CG/i && $bam !~ /20130415/) {  ## When a bas is not found in the same directory as BAM (usually ther archive staging area)
	my $bas_basename = basename($fo->name) . ".bas";
	my $bas_base = $fa->fetch_by_filename($bas_basename);
	
	if ( $bas_base && @$bas_base > 1  ) {
		write_log($fo, $loga, "FARM: multiple bas files found for $bam. One is " . $bas_base->[0]->name . " None is in the same directory as its BAM");
		throw("Multiple bas files found. One is " . $bas_base->[0]->name . " None is in the same directory as its BAM\n");
	}	
	else { # create a bas, store in db as new or replace the one exists in db

		eval{
			my $bas = ReseqTrack::Tools::Bas->new (
					 -reference => $reference,
                     -input_files => $bam,
                     -md5 => $fo->md5,
                     -working_dir => $working_dir,
                     -need_tags=> 0,
                     -in_parent => 1,
                    );
			$bas->run; ## this will create a bas file in the same directory as the bam file, when -in_parent is set to 1; otherwise specify -output_dir
		};
		if($@){
    		print "Cannot generate bas file for $bam\n";
    		print "$@\n";
    		throw("Cannot generate bas file for $bam\n$@");
  		}		

		my $bas_type = $fo->type;
		$bas_type =~ s/BAM/BAS/i;

		if (!$bas_base || @$bas_base == 0 ) { ## when there was no old bas object			
			$bas = store_file($bas_name, $bas_type, $host, "Created a bas file for BAM that doesn't come with a bas");
		}	
		else { ## when there is an old bas object (likely one that has been archived on the ftp site)
			$bas = update_file($bas_name, $bas_type, $host, "Update bas file", $bas_base);
		}				   
	}	
}

######################################################################################################
### Touch create BAI file for unmapped BMAs if they are not uploaded; load it into the database.   ###
######################################################################################################

my $bai_name = $bam . ".bai";
my $bai = $fa->fetch_by_name($bai_name);

if(!$bai) {

	my $bai_basename = basename($fo->name) . ".bai";
	my $bais = $fa->fetch_by_filename($bai_basename); # this is to query with "name like "%$bai_basename"
	
	if (@$bais >= 2) {
		throw("Have multiple files associated with ".$bai_basename);
	}
	elsif ( @$bais == 1) {
		write_log($fo, $loga, "FARM: bai file " . $bais->[0]->name . " is not in the same directory as its BAM file");
		throw("bai file " . $bais->[0]->name . " is not in the same directory as its BAM file\n");
	}	
	else {	
		if ($bam =~ /unmapped/ ) { ## Sanger not uploading .bai files for unmapped bams because they are all size zero with the same md5
			`touch $bai_name`;		
			my $bai_type = $fo->type;
			$bai_type =~ s/BAM/BAI/;
			$bai = store_file($bai_name, $bai_type, $host, "Created a bai file for unmapped bams");	
		}
		else {
			write_log($fo, $loga, "FARM: No bai file found for mapped BAM file, it is not even in the dropbox");
			throw("No bai file found for mapped BAM file $bam, it is not even in the dropbox\n");
		}
	}		
}	

SKIP:

######################################################################################################
### Check md5sum and archive files that have passed the check									   ###
######################################################################################################			
my $archive_list = $working_dir . "/" . basename($fo->name) . ".tmp_archive_list." . $time_stamp;

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

sub store_file {
	my ($f_path, $f_type, $host, $comment) = @_;
	my $hist_a = $db->get_HistoryAdaptor();
	my $f_objs = create_objects_from_path_list([$f_path], $f_type, $host); #$files is a reference to an array of file objects (in this case, only one element)
	my $f_obj = $f_objs->[0]; #to get the first and only file object
	my $md5 = run_md5($f_path);
	$f_obj->md5($md5); #set md5
	if ($run) {
		$fa->store($f_obj);	
		my $history = ReseqTrack::History->new(
			-other_id 	=> $f_obj->dbID,
			-table_name => 'file',
			-comment 	=> $comment, 
		);
		$hist_a->store($history);
		$f_obj->history($history);
	}
	else {
	    write_log($fo, $loga, "FARM: -run tag is not set so file created is not stored in database.");
		throw("-run tag is not set so file created is not stored in database.\n");
	}
	return $f_obj;
}	


sub update_file {
	my ($f_name, $f_type, $host, $input_comment, $f_base_objs) = @_;
	
	my $md5 = run_md5($f_name);	
	my $new_f_object = ReseqTrack::File->new
		(
			-adaptor => $fa,
			-dbID => $f_base_objs->[0]->dbID,
			-name => $f_name,
			-md5 => $md5,
			-host => $host,
			-type => $f_type,
		);
			
		my $history_ref = $f_base_objs->[0]->history;
		my $history;
			
		if (!$history_ref || @$history_ref == 0) {     	  
			$history = ReseqTrack::History->new(
				-other_id => $f_base_objs->[0]->dbID,
				-table_name => 'file',
				-comment => $input_comment, 
			);
			$new_f_object->history($history);
		}
		else {
			my $comment = calculate_comment($f_base_objs->[0], $new_f_object); 
				
			if (!$comment) {
				$comment = "fiddling files, no comment\n";
			}	
			$history = ReseqTrack::History->new(
				-other_id => $f_base_objs->[0]->dbID,
				-table_name => 'file',
				-comment => $comment,
			);
			$new_f_object->history($history);	
		}
			
		if ($run) {
			$fa->update($new_f_object, 1, 1); # the second 1 is allow change name		
		}
		else {
			write_log($fo, $loga, "FARM: -run tag is not set so bas file created is not stored in database.");
			throw("-run tag is not set so bas file created is not stored in database.\n");
		 }		
	return 1;
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
	-out				directory where output files will be written, default is current dir
	-working_dir		directory where temporary files will be written
	-reference			path to a reference genome fastq file which is used for bas creation
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

perl $ZHENG_RT/scripts/qc/bam_md5check_and_archive.pl -dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-bam /nfs/1000g-work/G1K/archive_staging/test/NA06985/alignment/NA06985.chrom10.ILLUMINA.bwa.CEU.low_coverage.2011b.bam \
-reference /nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/hs37d5/hs37d5.fa \
-working_dir /nfs/1000g-work/G1K/scratch/zheng/tmp

Older version of reference I used was /nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/bwa/grc37/human_g1k_v37.fa


