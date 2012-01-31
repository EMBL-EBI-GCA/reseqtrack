#!/sw/arch/bin/perl

use strict;

#use lib '/homes/zheng/reseq-personal/zheng/lib/reseqtrack_hzb/modules/'; #svn check out the ReseqTrack 
use ReseqTrack::DBSQL::RejectLogAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use File::Basename;
use File::Path;
use Getopt::Long;

$| = 1; 

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $help,
    $release_date,
    %status,
	%count_by_status,
   );

my $out_dir 				= "./";
my $myHost	 				= "1000genomes.ebi.ac.uk";
my $type 					= "BAM";
my $file_in_staging_flag 	= 0;

&GetOptions(
  'out:s'		=> \$out_dir,
  'dbhost=s'    => \$dbhost,
  'dbname=s'    => \$dbname,
  'dbuser=s'    => \$dbuser,
  'dbpass=s'    => \$dbpass,
  'dbport=s'    => \$dbport,
  'help!'		=> \$help,
  'host_name=s'	=> \$myHost,
  'date=s'		=> \$release_date,
);

if ($help) {
	help_info();
}

if (!$release_date && !$myHost) {
	throw("Please provide release date and host name; date has to be one of the seq index release dates in the formate of yyyymmdd eg 20101123\n");
}	

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host 	=> $dbhost,
  -user 	=> $dbuser,
  -port 	=> $dbport,
  -dbname 	=> $dbname,
  -pass 	=> $dbpass,
);

my $fa = $db->get_FileAdaptor;
#my $file_objs_in_this_release = $fa->fetch_by_type($type);  #for all sanger 20100901 release
my $file_objs_in_this_release = $fa->fetch_all_like_path("%$release_date%");

throw("No file object is found; please check the release date you entered, it needs to formatted as YYYYMMDD\n") if (!$file_objs_in_this_release || @$file_objs_in_this_release == 0);

my $loga = $db->get_RejectLogAdaptor;

foreach my $fo (@$file_objs_in_this_release) {
	
	next if ( $fo->type !~ /BAM|BAS|BAI/i );
	next if ($fo->type =~ /WITHDRAWN|TEST/i);  

	if ($myHost =~ /sanger/i) {
		next if ($fo->name !~ /bwa|ssaha|smalt/i );
		#next if ($fo->type =~ /EXOME/i);
	}
	elsif ($myHost =~ /tgen/i ) {
		next if ($fo->name !~ /SOLID/i);
		next if ($fo->type =~ /EXOME/i);
	}
	elsif ($myHost =~  /ncbi/i ) {
		next if ($fo->name !~ /mosaik/i);
		next if ($fo->type =~ /EXOME/i);
	}		
	elsif ($myHost =~  /broad/i ) {
		next if ($fo->name !~ /bwa/i);
		next if ($fo->type !~ /EXOME/i);
	}	
	elsif ($myHost =~  /baylor/i ) {
		next if ($fo->name !~ /SOLID/i);
		next if ($fo->type !~ /EXOME/i);
	}	
	elsif ($myHost =~  /boston_college/i ) {
		next if ($fo->name !~ /mosaik/i);
		next if ($fo->type !~ /EXOME/i);
	}		
	
	my $ha 			= $db->get_HostAdaptor;
	my $host_obj 	= $ha->fetch_by_name($myHost);
	my $dropbox_dir = $host_obj->dropbox_dir;
	my $full_path 	= $fo->name;
	my $basename 	= basename($full_path);
	
	my $log_objs 	= $loga->fetch_by_file_id($fo->dbID);
	
	if (!$log_objs || @$log_objs == 0) {
	#	warning("No reject_log is found for file " . $fo->name . "\n");
	}
	elsif (@$log_objs >1) {
		throw("File " . $fo->name . "has more than one log\n");
	}
	
	my $log_obj = $log_objs->[0];
	
	if ($fo->name =~ /\/nfs\/1000g-work\/G1K\/drop\/reject/) {
		$status{$fo->name} = "Transferred:failed_qa\t" . $log_obj->reject_reason;
		$count_by_status{"Transferred:failed_qa"}++;
		#### capture md5 check failures, bas content inconsistencies, files with wrong sizes etc.
	}
	elsif ($fo->name =~ /\/nfs\/1000g-work\/G1K\/archive_staging\// ) {
		if ($log_obj) {
			$status{$fo->name} = "Transferred:in_staging\t" . $log_obj->reject_reason;
		}
		else {
			$status{$fo->name} = "Transferred:in_staging";
		}	
		$count_by_status{"Transferred:in_staging"}++;
		$file_in_staging_flag = 1;
		#### failures due to missing bas or bai files can be found this way; the farm job failed but bam is not put in the reject bin
	}	
	elsif ($fo->name =~ /\/vol1\/ftp\//) {
		$status{$fo->name} = "Transferred:done\t";
		$count_by_status{"Transferred:done"}++;
	}
	elsif ($fo->host_id != 1) { #if the host name is remote, not  "1000genome.ebi.ac.uk" 
		$status{$fo->name} = "Loaded:in_db\t";
		$count_by_status{"Loaded:in_db"}++;		
=head
		my $derived_file_path_in_dropbox = `find $dropbox_dir -name $basename -print`; 
		chomp $derived_file_path_in_dropbox;
		
		if (-e $derived_file_path_in_dropbox)  {
			$status{$fo->name} = "Transferred:in_process\t";
			$count_by_status{"Transferred:in_process"}++;
		}
		else {
			my $file_timestamp = $fo->updated; # 2009-11-23 09:47:55
			my @bits = split (/\s+/, $file_timestamp);
			my ($yr, $mon, $day) = split (/-/, $bits[0]);
			my ($hr, $min, $sec) = split (/:/, $bits[1]);
				
			my $days_elapsed = ( time() - timelocal($sec, $min, $hr, $day, $mon-1, $yr-1900 ) )  / 86400; 
				
			if ( $days_elapsed > 7) {  
			##### If the file has been in the db for more than 7 days.  FIXME magic number, need reality check
				$status{$fo->name} = "Transfer never happened\t";
				$count_by_status{"Transfer never happened"}++;
			}	
			else {
				$status{$fo->name} = "To be transferred\t";
				$count_by_status{"To be transferred"}++;
			}
		}
=cut
	}
	else {
		$status{$fo->name} = "Status not known\t";
		$count_by_status{"Status not known"}++;
	}			
}

print "SUMMERY for $myHost center, $release_date release\n";
foreach my $possible_status (keys %count_by_status) {
	print "$possible_status\t$count_by_status{$possible_status}\n";
}

#### if there are files in the staging area, find out the farm job status. ####
#### While some failed farm jobs (for instance failed due to md5 check) would result in files being put in the reject bin #### 
#### other failed jobs (caused by out-of-memory and other reasons) do not result in having the files moved into the reject bin ####

if ( $file_in_staging_flag == 1 ) { 
#if ($file_in_staging_flag == 0) {	
	#### FIXME: use ==1 in real cases
	my $ja = $db->get_JobAdaptor;
	my $ca = $db->get_EventCompleteAdaptor;
	my $aa = $db->get_EventAdaptor;
	
	my ($jobs, 
		$analyses, 
		$input_hash, 	
		%status_count,
		%logic_status_count,
		%input_string,
		%event_hash
	);
	
	$jobs       = $ja->fetch_all;
	$analyses   = $aa->fetch_all;
	
	my $completed_event = get_completed_event_hash($db);

	foreach my $analysis (@$analyses) {
	    $event_hash{$analysis->name} = $analysis;
	}
	
	foreach my $job (@$jobs) {
	    if (!$status_count{$job->current_status}) {
	        $status_count{$job->current_status} = 1;
	    } else {
	        $status_count{$job->current_status}++;
	    }
	    if (!$logic_status_count{$job->event->name}) {
	        $logic_status_count{$job->event->name} = {};
	    }
	    if (!$logic_status_count{$job->event->name}{$job->current_status}) {
	        $logic_status_count{$job->event->name}{$job->current_status} = 1;
	    } else {
	        $logic_status_count{$job->event->name}{$job->current_status}++;
	    }
	}
	
	print "-----------\nBam QA farm jobs status - on-going\n";
	foreach my $name (keys(%logic_status_count)) {
		next if ($name !~ /bam/i ); # so won't print events that have nothing to do with BAM releases
        foreach my $status (keys(%{$logic_status_count{$name}})) {
            print $name. " "
              . $status . " "
              . $logic_status_count{$name}->{$status} . "\n----";
              
            #if ($status eq "FAILED") {  
            #	print "\nUse the following sql statement to find which jobs have failed:\n";
            #	print "**** select job.stderr_file, job.exec_host from job, job_status where job.job_id = job_status.job_id and is_current = 'y' and status = 'FAILED' and event_id =23****\n\n";  
        	#}
        }
    }
    
    print "\nBam QA farm jobs - completed\n";
    foreach my $logic (keys(%$completed_event)) {
        next if ($logic !~ /bam/i);
        my $event = $event_hash{$logic};
        my $count = $completed_event->{$logic}->{$event->type};
        print $logic. "\t" . $event->type . "\t" . $count . "\n";
    }
    print "-----------\n";
}	

print "\nDETAILS\n";
print "File\tStatus\tReject_reason\n";	
foreach my $file ( keys %status) {
	print "$file\t$status{$file}\n";
}

####### SUBS ########	
sub get_completed_event_hash {
    my ($db) = @_;
    my $sql =
"select event.name, event.type, count(distinct(other_id)) from event, event_complete where event.event_id = event_complete.event_id group by event.name , event.type";
    my $sth = $db->dbc->prepare($sql);
    $sth->execute;
    my %hash;
    while (my ($name, $type, $count) = $sth->fetchrow) {
        $hash{$name}->{$type} = $count;
    }
    return \%hash;
}

=head

perl ~/reseq-personal/zheng/bin/check_bam_release_status.pl -dbhost mysql-g1kdcc -dbname g1k_archive_staging_track -dbuser xxxx -dbpass xxxx -dbport 4197 -host_name sanger -date 20101123 > status

TEST:
perl ~/reseq-personal/zheng/bin/check_bam_release_status.pl -dbhost mysql-g1kdcc -dbname zheng_automation_test -dbuser xxxx -dbpass xxxx -dbport 4197 -host_name sanger -date 2011b
