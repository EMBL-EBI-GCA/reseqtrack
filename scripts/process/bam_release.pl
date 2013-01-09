#!/sw/arch/bin/perl -w

use strict;

use ReseqTrack::DBSQL::DBAdaptor; 
use ReseqTrack::DBSQL::RejectLogAdaptor;
use ReseqTrack::DBSQL::HostAdaptor;  
use ReseqTrack::Tools::BamUtils;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use File::Basename;
use File::Path;
use Getopt::Long;
use Time::Local;

#use lib '/nfs/1000g-work/G1K/work/zheng/reseq-personal/zheng/lib';
#use myTIME;
use ReseqTrack::Tools::myTIME;

$| = 1; 

my $start_run = time();

my ($time_stamp, $month_stamp, $day_stamp) = ReseqTrack::Tools::myTIME::get_time();

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $help,
    $run_dir,
    $run,
    $verbose,
    $analysis_grp,
    $desired_host,
   );

my $move_to_dir = "/nfs/1000g-work/G1K/archive_staging/ftp/data";
my $output_dir = ".";
my $run_farm_job_flag = 0; #when all files have been moved to archive_staging area; just need to kick off farm job 

&GetOptions(
  'dbhost=s'     		=> \$dbhost,
  'dbname=s'     		=> \$dbname,
  'dbuser=s'     		=> \$dbuser,
  'dbpass=s'     		=> \$dbpass,
  'dbport=s'     		=> \$dbport,
  'help!'		 		=> \$help,
  'move_to_dir=s'		=> \$move_to_dir,
  'run!'				=> \$run,
  'verbose!'			=> \$verbose,
  'out=s'				=> \$output_dir,
  'kick_off_farm_job:s'	=> \$run_farm_job_flag, #value is ncbi, sanger, tgen, broad, baylor, boston_college
  'analysis_grp:s'		=> \$analysis_grp, #value is low_coverage or exome;  this is to provide analysis group when -kick_off_farm_job is used 
  'host:s'				=> \$desired_host, #this is to only run a specific host
);

if ($help) {
	help_info();
}

if ($run_farm_job_flag) {
	if ( $run_farm_job_flag !~ /ncbi|sanger|tgen|baylor|boston_college|broad/ ) {
		throw("-kick_off_farm_job option has to use one of the follows as input: ncbi, sanger, tgen, baylor, boston_college, broad\n"); 
	}
}	

print "\n*****ARE YOU RUNNING IN A PRODUCTION SERVER SUCH AS EBI-002 or PG-TRACE-001?********\n\n";

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );
    
my $ha = $db->get_HostAdaptor;
my $remote_hosts = $ha->fetch_all_remote();

throw("No remote host found, please check which database you are using") if (!$remote_hosts);

my $loga = $db->get_RejectLogAdaptor;
my $fa = $db->get_FileAdaptor;

$db->dbc->disconnect_when_inactive(1);

my $db_parameters = "-dbhost $dbhost -dbuser $dbuser -dbpass $dbpass -dbport $dbport -dbname $dbname";

$output_dir =~ s/$\///g;
	 
foreach my $host ( @$remote_hosts ) {
	
	if ( $run_farm_job_flag ) {
		next if ($run_farm_job_flag ne $host->name);
	}
	
	#print "host is " . $host->name . ", desired host is $desired_host\n";
	if ($desired_host) {
		next if ($host->name ne $desired_host);
	}
		
	my %dir_hash;	
	my %dir_file_hash;

	my $in_process_flag = 0;	# this is to tell whether all files in db and all files in dropbox have been processed; when all done, $flag should remain 0   
	my $fail_flag 		= 0;	# if any file stuck during transfer after 2 days, this flag will be set to 1
	my $move_flag 		= 0;	# if any files have been processed successfully

	my $host_id 		= $host->dbID;
	my $host_name 		= $host->name;
	my $dropbox_dir 	= $host->dropbox_dir;
	
	my $log_file =  $output_dir . "/log." . $host_name . "." .$time_stamp;
	open (my $log, ">", $log_file) || throw("Cannot open log file $log_file\n");
	
	my $withdraw_file = $output_dir . "/withdraw_list." . $host_name . "." . $time_stamp;
	open (my $withdrawn_list_fh, ">", "$withdraw_file") || throw ("Cannot open $withdraw_file\n"); 
	
	#### To make sure all bam/bai/bas files in the database make a corresponding entry in the dropbox 
	#### To make sure that file size in dropbox is the same as the one in the db  
	#### To move a qualified file to a specified dir and later insert it into a collection
	  
	my $files = $fa->fetch_by_host($host_id);
	#warning ("No file object is found; please check the host $host_name\n") if (!$files || @$files == 0 );
	
	my %analysis_group;
	
	if ($analysis_grp) {
		$analysis_group{$analysis_grp} = 1; # this is to memorize the user input analysis group
	}
	
	foreach my $file ( @$files ) {
		print "Host $host_name: ". $file->name . "\n" if ($verbose);
		my $original_file_name = $file->name;
		
		my $size_in_db = $file->size;
		my $basename = basename($original_file_name);
		my $file_name = `find $dropbox_dir -name $basename -print`; 
		#print "find file: $file_name\n";
		chomp $file_name;
		
		if ($file_name =~ /\n/) { # if after chomp, still new line remains, that means >1 files are pulled out
			write_log($file, $loga, "PREP: more than one file have the same name");
			print STDERR "PREP: more than one file have the same name $file_name\n";
			next;
		}	

		if ( $file->type =~ /EXOME/ ) {
			$analysis_group{'exome'} = 1;
		}
		elsif ( $file->type =~ /CG/ || $file->type =~ /HIGH_COV/ ) {
			$analysis_group{'high_coverage'} = 1;
		}	
		else {
			$analysis_group{'low_coverage'} = 1;
		}
				
		my $aspx = $file_name . ".aspx" if ($file_name);
		
#		if ($file->type =~ /BAM/ || $file->type =~ /BAI/ || $file->type =~ /BAS/ ) {
		if ($file->type =~ /BAM/ || $file->type =~ /BAI/ || $file->type =~ /BAS/ || $file->type =~ /CG_/ ) {
			if ( $file_name && -e $file_name) {
				my $size = -s $file_name;
				if ( $size != $size_in_db && ( ( -M $file_name ) > 99) ) {	# if the file size is different from what is the db 
																			# and the file has been 100 days or more old 
					$fail_flag = 1;
					write_log($file, $loga, "PRE-PROCESSING: uploading failed? file in dropbox has different size as the one in db (dropbox $size, db $size_in_db) after 100 days");
					move_bam_to_trash($db, $file, $file_name, $run);
				}
				elsif ( $size == 0 && $size_in_db == 0) {
					$fail_flag = 1;
					write_log($file, $loga, "PRE-PROCESSING: uploaded file with size 0"); 
				}
				elsif (	(-s $file_name) != $size_in_db ) {					# if the file size is different from what is the db 
					$in_process_flag = 1;
					write_log($file, $loga, "WARNING: upload on going? file in dropbox has different size as the one in db (dropbox $size, db $size_in_db)");
				}
				#elsif ( $aspx && (-e $aspx) &&  ( (-M $aspx ) > 2) ) {
				#	$fail_flag = 1;
				#	write_log($file, $loga, "PREP: uploading failed? file in dropbox still has an aspx file after 3 days");
				#	move_bam_to_trash($db, $file, $file_name, $run);
				#}		
				#elsif ($aspx && (-e $aspx) ) {
				#	$in_process_flag = 1;
				#	write_log($file, $loga, "WARNING: uploading on-going? file in dropbox still has an aspx file");
				#}		
				else { 

					if ( $file->type =~ /BAS/ && check_bas($file_name) == 1) {
						write_log($file, $loga, "PREP: bas file has content inconsistency");
						$fail_flag = 1;
						move_bam_to_trash($db, $file, $file_name, $run);
					}	
					else {	
						my ($new_directory, $new_f_obj) = check_name_and_move_file($file, $file_name, $move_to_dir, $host, $run); 	
						# $new_directory is where the file is after loaded
						if ($new_directory && $run) {
							$move_flag = 1;
							$dir_hash{$new_directory} = 1; 
							push @{$dir_file_hash{$new_directory}}, $new_f_obj;
						}		
						write_log($file, $loga);

=head

						if ($run) {
							if ( $new_f_obj->type =~ /BAS/ && check_bas($new_f_obj->name) == 1) {
								write_log($file, $loga, "PREP: bas file has content inconsistency");
								$fail_flag = 1;
								move_bam_to_trash($db, $new_f_obj, $new_f_obj->name, $run);
							}
						}		

=cut					

						####### FIXME: use above lines when doing changing file names for TGEN	
					}
				}
			}	
			else { # when a file is not found in the dropbox
				my $file_timestamp = $file->updated; # 2009-11-23 09:47:55
				my @bits = split (/\s+/, $file_timestamp);
				my ($yr, $mon, $day) = split (/-/, $bits[0]);
				my ($hr, $min, $sec) = split (/:/, $bits[1]);
				
				my $days_elapsed = ( time() - timelocal($sec, $min, $hr, $day, $mon-1, $yr-1900 ) )  / 86400; 
				
				if ( $days_elapsed < 15) {  ##### FIXME, 8 is a magic number, need reality check
					$in_process_flag = 1;
					write_log($file, $loga, "WARNING: File does not exist in dropbox, it might not have been loaded yet!");
				}
				else {
					$fail_flag = 1;
					write_log($file, $loga, "PREP: File does not exist in dropbox after 7 days it was loaded in db, transfer failed?");
				}		
			}
		} # end of checking type	
		else {
			print STDERR "File $original_file_name is not a BAM/BAS/BAI, ignore!\n";
			write_log($file, $loga, "File is not of type BAM/BAS/BAI, ignore!");
		}			
	} # end of foreach file	
	
	#### To make sure all files in the dropbox is in the db as well ####

	$dropbox_dir =~ s/\/$//; # remove the last back slash
	my ($files_in_dropbox, $hash) = list_files_in_dir($dropbox_dir); # hash key is dir name, value is an array of file names
	
	if (!$files_in_dropbox) {
		print STDERR "WARNING: No file found in dropbox $dropbox_dir, perhaps they have all been moved to $move_to_dir\n";
	}	
	
	foreach my $dr (keys %$hash) {
		my %h = %$hash;
		#print "my dir is $dr\n";
		#print "files in dropbox $host_name are " . join("\n", @{$h{$dr}}) . "\n";
		foreach my $f (@{$h{$dr}}) {
			next if ($f =~ /aspx/);   
			next unless ($f =~ /bam/);   
			next if ($f =~ /md5/);
			my $full_file_path = $dr . "/" . $f;
			#my $fo = $fa->fetch_by_name($full_file_path);	
			my $fos = $fa->fetch_by_filename($f);
			if (!$fos || @$fos == 0 ) {
				$fail_flag = 1;
				print $log "WARNING: BAM file $full_file_path in dropbox is not loaded in database.\n";
			}	
			else {
				if ( $fos->[0]->name =~ /vol1/ ||  $fos->[0]->name =~ /staging/ ) {
					print $log "WARNING: BAM file $full_file_path in dropbox already exists on staging area or on ftp site; is it an updated version?\n";
				}
				else {
					print "in process: " . $fos->[0]->name . "\n";
				}
			}
		}	
	}
	
	if ($move_flag == 1) {			
		foreach my $alignment_dir (keys %dir_hash) {
			create_and_load_collection($run, $db, $alignment_dir, $log, $withdrawn_list_fh, $withdraw_file);
		}
	}

	my @groups = keys %analysis_group;
	
	foreach $analysis_grp ( @groups ) {
	
		my $command = "perl /nfs/1000g-work/G1K/work/zheng/reseqtrack/scripts/event/run_event_pipeline.pl ";
		$command .= "-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 ";
#		$command .= "-dbhost mysql-g1kdcc-public -dbname zheng_var_call -dbuser g1krw -dbpass thousandgenomes -dbport 4197 "; #### FIXME, change to the above line after testing
		$command .= "-once -runner /nfs/1000g-work/G1K/work/zheng/reseqtrack/scripts/event/runner.pl "; 	
		
		## FIXME, tgen will provide exome bams
		
		if ( ($host_name eq "sanger" || $host_name eq "tgen") && defined $analysis_grp && $analysis_grp eq "low_coverage" ) {
			$command .= "-name bam_release & ";
		}	
		elsif ( $host_name eq "sanger" && defined $analysis_grp && $analysis_grp eq "exome" ) {
			$command .= "-name exome_bam_release & ";
		}	
		elsif ( $host_name eq "sanger" && defined $analysis_grp && $analysis_grp eq "high_coverage" ) { ## This is a hack to release the high cov CEU trio BAMs Sanger re-mapped
			$command .= "-name ncbi_bam_release & ";
		}	
		elsif ( $host_name eq "tgen" && defined $analysis_grp && $analysis_grp eq "exome" ) {
			$command .= "-name exome_bam_release & ";
		}
		#elsif ($host_name eq "baylor" && defined $analysis_grp && $analysis_grp eq "exome" ) {
		#	$command .= "-name exome_bam_release & ";
		#}	
		elsif ($host_name eq "ncbi") {
			$command .= "-name ncbi_bam_release & ";
		}	
		elsif ($host_name eq "boston_college" && defined $analysis_grp && $analysis_grp eq "exome" ) {
			$command .= "-name exome_bam_release_boston_college & ";
		}	
		#elsif ( $host_name eq "broad") { 
		#	$command .= "-name exome_bam_release_broad & ";
		#}	
		elsif ( ! $analysis_grp  && $run_farm_job_flag ) {
			warning("For host $host_name no analysis group has been parsed out, please provide -analysis_grp to -kick_off_farm_job\n");
		}	
		else {
			warning("Problem - don't know what to do with host $host_name and analysis group $analysis_grp"); 
		}	
		
		#print "farm command is $command\n";
		
		if ( $move_flag == 1  || $run_farm_job_flag eq $host_name ) {	# sometimes the files have been moved to archive staging but 
																		# farm jobs did not get sent due to various reasons
																		# the -kick_off_farm_jobs option allows run_event_pipeline to run
			
			#### after collection is created for a given file, kick off farm jobs to run check md5s and archive files. 
			#### for event run_bam_md5check_and_archive, any files that are not in archive_staging area will fail the run
			#### any bams whose bas and/or bai are not on the archive staging area will fail the run (except the NCBI bams)
			#### inside run_bam_md5check_and_archive, move failed files to reject bin
			
			print "Submitting farm jobs.........\n";												
			`$command`;

=head										
	
			foreach my $directory ( keys %dir_file_hash) {
				foreach my $f ( @{$dir_file_hash{$directory}} ) {
					my $file_to_archive = $f->name;
					next if ($file_to_archive =~ /bai/ || $file_to_archive =~ /bas/);
					#print "Input file name for the qa program is: $file_to_archive\n" if ($verbose);
					`/usr/bin/perl /homes/zheng/reseq-personal/zheng/bin/bam_md5check_and_archive.auto.pl $db_parameters -bam $file_to_archive -out $output_dir -verbose`;
					# FIXME: for test, did not use the -run option since zheng-automation_test is not capable of archiving.  
					# use the -run option after test!!  
				}
			}			
	
=cut
	###FIXME: comment out above lines after testing
		}
		else {
			print "No file is moved and processed this time for host $host_name, did you set the -run tag?\n";
		}
	}# end of foreach @groups
	
	
	if ($fail_flag == 1) {
		print "Some files failed transfer or some files in dropbox do not exist in database. Please check the log file $log_file and reject_log table in the db\n";
	}		
	elsif ($in_process_flag == 1 ) {
		print "Some files are still being transferred, cron job to sleep\n";
	}
	elsif ($run) {
		print "ALL files have been processed, can proceed to withdraw files and dump alignment index\n";
	}	
	else {
		print "No files have been processed. The reason can be no file has been uploaded or the -run flag is not set!\n";
	}
	  
	close $withdrawn_list_fh;
	close $log;
	`rm $log_file` if ( -s $log_file == 0);
	`rm $withdraw_file` if ( -s $withdraw_file == 0); 
	print "****** Finished processing host $host_name ******\n\n";
} #end of one host

####################
###### SUBS ########
####################

#### Create collections for loaded BAM files if -load_collection is set
sub create_and_load_collection {
	
	my ($run1, $db1, $path, $log, $withdrawn_list_fh, $withdraw_file) = @_;

	my $bam_cnt = 0;
	my $file_type;
	my $chrom;
	my %collection_hash_array; 	#key is the name of the collection, value is an array of file objs that belong to this collection
	my %collection_hash_hash;
	my %collection_chr; 		#key is collection name, second key is chromosomes that are being loaded	
	my %collection_date;
	
	$path =~ s/\/$//; 			#remove last slash 
	my $fuzzy_path = $path . "/%.bam";
	print $log "path is $fuzzy_path\n" if ($verbose);
			
	my $fa = $db1->get_FileAdaptor;
	my $ec = $db1->get_EventCompleteAdaptor;
	
	my $f_objs =  $fa->fetch_all_like_path($fuzzy_path); # this returns bai and bas objects as well!
	warning("No files are found in path $fuzzy_path\n") if (!$f_objs || @$f_objs==0);

	#my $f_objs  = $fa->fetch_all_like_path("/nfs/1000g-archive/vol1/ftp/pilot_data/data/%/alignment/NA07347%.bam");
	
	foreach my $fo (@{$f_objs}) {
		my $file_path = $fo->name;
		if ($fo->type =~ /BAM/i) {
			$file_type = $fo->type;			
			my ($collection_name, $sample, $platform, $algorithm, $project, $analysis, $chrom, $date) = get_collection_name_from_file_name($file_path);
				
			push @{$collection_hash_array{$collection_name}}, $fo;
			$collection_chr{$collection_name}{$chrom} =1;
			$collection_date{$collection_name} = $date;
			$collection_hash_hash{$collection_name}{$fo->filename} = 1;
			$bam_cnt++;
		}	
	}
	
	print $log "number of bam files is $bam_cnt for one individual\n" if ($verbose);
	goto SKIP if ($bam_cnt == 0);
		
	foreach my $coll_name ( keys %collection_hash_array) {	
		
		my $ca = $db->get_CollectionAdaptor;
		my $exist_collection = $ca->fetch_by_name_and_type($coll_name, $file_type);
		
		my $fo_to_be_considered_to_withdraw;
		my @bams_to_remove;
		 
		if ($exist_collection) {
			
			###### when update happens, new BAM files with new date are loaded, a withdraw file list 
			###### will report old files that should be withdrawn (basically old files with the same collection name with the new load

			$fo_to_be_considered_to_withdraw = $exist_collection->others; 
			foreach my $old_bam (@$fo_to_be_considered_to_withdraw) {
				
				my $old_bam_name = $old_bam->name;
				my $old_bam_basename = $old_bam->filename;
				
				#if ( $old_bam_name =~ /technical\/working/ && $old_bam->type ne "EXOME_BI_BAM") { # no need to withdraw files in tech/working dir, 
																								   #just remove them from collection_group table
				if ( $old_bam_name =~ /technical\/working/ ) {
					push @bams_to_remove, $old_bam;
					next;
				}
					
				my ($sample2, $platform2, $algorithm2, $project2, $analysis2, $chrom2, $date2) = CHECK_AND_PARSE_FILE_NAME($old_bam_name);
				
				if ($old_bam_name =~ /\/nfs\/1000g-archive\/vol1\// &&
					$collection_chr{$coll_name}{$chrom2}  &&									
				    $date2 != $collection_date{$coll_name} ) {

				#### only withdraw files that have been loaded on the ftp site; the ones that are in staging area might be there as part of an unfinished load
				#### only withdraw files that have a different date
				#### only withdraw files for chromosomes that are being loaded, not a blanket withdraw							
					
					#print $log "Inside if loop: new date: $date2, old date $collection_date{$coll_name}\n";			
					push @bams_to_remove, $old_bam;
		
					my $old_bas_name = $old_bam_name . ".bas";
					my $old_bas_obj = $fa->fetch_by_name($old_bas_name);
					
					my $old_bai_name = $old_bam_name . ".bai";
			
					my @old_bam_name_bits = split (/ftp/, $old_bam_name);
					my $withdraw_bam_name = $old_bam_name_bits[0] . "withdrawn" . $old_bam_name_bits[1];
					my $withdraw_bas_name = $withdraw_bam_name . ".bas";
					my $withdraw_bai_name = $withdraw_bam_name . ".bai";
					
					print $withdrawn_list_fh "$old_bam_name\t$withdraw_bam_name\n";  
					print $withdrawn_list_fh "$old_bas_name\t$withdraw_bas_name\n" if ($old_bas_obj); ## NCBI release does not have bas files for now;
					print $withdrawn_list_fh "$old_bai_name\t$withdraw_bai_name\n";
				}
				
				## When a file is replaced with another file with the same name, need to remove the 
				## file id from the event_complete table so the md5 check can be run on the farm on the new file
				if ($collection_hash_hash{$coll_name}{$old_bam->filename}) { #if the old file in the collection have the same name as the newly uploaded one
					my $old_bam_id = $old_bam->dbID;					
					my $old_bam_records_in_ec = $ec->fetch_by_other_id($old_bam_id, "file"); #file is the table name
					if (!$old_bam_records_in_ec || @$old_bam_records_in_ec == 0) {
						print "no event_complete record is found for " . $old_bam->filename ."\n" if ($verbose);
					}
					else {
						foreach my $ec_record (@$old_bam_records_in_ec) {
							$ec->remove($ec_record);
							print "event_complete record removed for replacement file of the same name " . $old_bam->filename ."\n" if ($verbose);
						}
					}		
				}					
			}#end of foreach old bam
		}

		my $collection =  ReseqTrack::Collection->new(
		  -name => $coll_name,
		  -others => \@{$collection_hash_array{$coll_name}},
		  -type => $file_type,  
		  -table_name => 'file',
		);
	
		if($run1 ) {
			$ca->store($collection); ## the store function checks for exisiting collection and update them
			my $new_collection = $ca->fetch_by_name_and_type($coll_name, $file_type);			
			$ca->remove_others($new_collection, \@bams_to_remove) if ($exist_collection && @bams_to_remove > 0);
			
			print $log "Collection $coll_name is loaded\n" if ($verbose && @bams_to_remove == 0);
			print $log "Collection $coll_name is updated and the old files on ftp site are listed in $withdraw_file\n" if ($verbose && @bams_to_remove > 0);
		
		}
		elsif (!$run1) {
			print $log "Collection is not loaded. Please set -load_collection if you like to load them\n" if ($verbose);
		}
	}
	SKIP:
	return ($bam_cnt);
}	

			
sub check_name_and_move_file {
	
	my ($file, $full_name, $move_to_dir, $host, $run) = @_; #$full_name is actual file path
		
	my $filen = basename($full_name);
	my @tmp = split(/\./, $filen);
	my @tmp2 = split(/\_/, $tmp[0]);
	my @tmp3 = split(/-/, $tmp2[0]);
	my $ind = $tmp3[0];	
	
	$move_to_dir =~ s/\/$//;
	
	if ( $file->type =~ /BAM/i) {
		CHECK_AND_PARSE_FILE_NAME($filen);
	}		
	
	my $new_dir;
	
	#### FIXME: make sure these rules and file type assignment rules still apply to phase 2 data!!
	if ($filen =~ /mosaik/ && $filen !~ /exome/ ) { # For NCBI bams
		$move_to_dir = "/nfs/1000g-work/G1K/archive_staging/ftp/technical/ncbi_varpipe_data";
		$new_dir = $move_to_dir . "/alignment/" . $ind . "/";
	}
	elsif ( $filen =~ /exome/i ) {
		
		#if ( ($filen =~ /bwa/i && $host->name =~ /sanger/i) || ( $filen =~ /bfast/i && $host->name =~ /baylor/i) )  { # FOR sanger exome and Baylor exome
		#	 $new_dir = $move_to_dir . "/" . $ind . "/exome_alignment/";
		#}
		if ( ($filen =~ /bwa/i && $host->name =~ /sanger/i) || ( $filen =~ /bfast/i && $host->name =~ /tgen/i) )  { # FOR sanger exome and tgen exome
			 $new_dir = $move_to_dir . "/" . $ind . "/exome_alignment/";
		}		
		elsif ($filen =~ /mosaik/i) { # FOR Boston College exome 
			$move_to_dir = "/nfs/1000g-work/G1K/archive_staging/ftp/technical/other_exome_alignments";
			$new_dir = $move_to_dir . "/" . $ind . "/exome_alignment/";
		}
		else {
			throw("Cannot tell to where to move file $filen on host " . $host->name);
		}	
				
#		if ($filen =~ /bwa/i || $filen =~ /solid\.mosaik/i)  { # For Broad EXOME bwa illumina BAMs or Boston college SOLID mosaik BAMs
#			$move_to_dir = "/nfs/1000g-work/G1K/archive_staging/ftp/technical/other_exome_alignments";
#			$new_dir = $move_to_dir . "/" . $ind . "/exome_alignment/";
#		}
#		else { 
#			$new_dir = $move_to_dir . "/" . $ind . "/exome_alignment/";
#		}

	}	
	elsif ( $filen =~ /COMPLETE_GENOMICS/ ) {
		$new_dir = $move_to_dir . "/" . $ind . "/cg_data/";
	}	
	else {
		$new_dir = $move_to_dir . "/" . $ind . "/alignment/";
	}		
		
	unless (-e $new_dir) {
		mkpath($new_dir);
	}
			
	my $new_file_path = $new_dir . $filen;
				
	`mv $full_name $new_file_path` if ($run);	 	
	my $exit = $?>>8;
	throw("mv failed\n") if ($exit >=1);
		
	print "old file path is $full_name and new file path is $new_file_path\n" if ($verbose);		
	
	my $new_host = $ha->fetch_by_name("1000genomes.ebi.ac.uk");
					
	my $new_file_object = ReseqTrack::File->new
 	   (
	      -adaptor => $fa,
	      -dbID => $file->dbID,
	      -name => $new_file_path,
	      -md5 => $file->md5,
	      -host => $new_host,
	      -type => $file->type,
	      -size => $file->size,
	      -created => $file ->created
	   );
			     
	my $history_ref = $file->history;
	my $history;
			 
	if (!$history_ref || @$history_ref == 0) {     	  
		$history = ReseqTrack::History->new(
			-other_id => $file->dbID,
			-table_name => 'file',
			-comment => "Rename or move file to a different location while keep md5 unchanged", 
		);
		$new_file_object->history($history);
	}
	else {
		my $comment = calculate_comment($new_file_object, $file); 
		if (!$comment) {
			$comment = "fiddling data, no comment\n";
		}	
		$history = ReseqTrack::History->new(
			-other_id => $file->dbID,
			-table_name => 'file',
			-comment => $comment,
		);
		$new_file_object->history($history);
	}
	
	$fa->update($new_file_object, 1, 1) if($run); # the second 1 is allow change name
	            	
	return ($new_dir, $new_file_object);
}
			    
sub help_info {

 exec('perldoc', $0);

}


#############################################################################################

=pod

=head1 NAME

 ~/ReseqTrack/scripts/process/bam_release.pl

=head2 Required arguments:

	-dbhost, 			the name of the mysql-host
	-dbname, 			the name of the mysql database
	-dbuser, 			the name of the mysql user
	-dbpass, 			the database password if appropriate
	-dbport, 			the port the mysql instance is running on, this defaults to 4197 the standard
            				 port for mysql-g1kdcc.ebi.ac.uk
		
	The script needs to be run on a production node such as ebi-002.

=head2 Optional arguments:

	-run				when this tag is used, files will be moved and farm jobs will be submitted to check and archive BAMs 
	-out				directory where the log files will be written, default is the run dir
	-move_to_dir		directory where the files will be moved to before archiving; 
						default is /nfs/1000g-work/G1K/archive_staging/ftp/data and other paths derived based on file types and names
	-kick_off_farm_job	value for this tag is a host name such as "sanger". 
						It will force the system to look for files that haven't been processed by farm job bam_release or 
						exome_bam_release .... Tag -analysis_grp is required when this option is used.
	-analysis_grp		This is needed when kick_off_farm_job is set, to indicate "exome", "low_coverage"					
	-verbose			default is off, set flag -verbose to print run logs
	-help				this makes the script print out its options
	

=head1 SYNOPSIS

 Once BAM and BAI file paths and md5s are loaded in database with appropriate host name, this script will look for files with remote host names and do the following:
 	- if files have been full transferred into the dropbox, then they will be moved to proper location in archive staging area (defined by -move_to_dir)
 	- a collection is created for each file or a group of files
 	- a withdraw list is produced if there are old files belong to the same collection 
 	- kick off a farm job that does the followings:
 		- if no bas file, create bas file and store in db
	 	- for unmapped BAM, if no bai file exist, create one for 
	 	- check md5 for each file
	 	- archive files that have passed md5 check

=head1 Example:

 perl ~/ReseqTrack/scripts/process/bam_release.pl -dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass xxxxxxxxxxxx -dbport 4197 -verbose -out /nfs/nobackup/resequencing_informatics/zheng/bam_release -run &
 OR
 perl /nfs/1000g-work/G1K/work/zheng/reseqtrack/scripts/process/bam_release.pl -dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass xxxxxxxxx -dbport 4197 -move_to_dir /nfs/1000g-work/G1K/archive_staging/test -verbose -out /nfs/1000g-work/G1K/scratch/zheng/exome_bam_release -kick_off_farm_job baylor -analysis_grp exome 
 
 TEST:
 perl ~/ReseqTrack/scripts/process/bam_release.pl -dbhost mysql-g1kdcc-public -dbname zheng_automation_test -dbuser g1krw -dbpass xxxxxxxxxx -dbport 4197 -move_to_dir /nfs/1000g-work/G1K/archive_staging/test -verbose -run
