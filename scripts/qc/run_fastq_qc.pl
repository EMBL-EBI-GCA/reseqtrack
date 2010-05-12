#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::StatisticsUtils;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use File::Basename;
use File::Path;
use File::Copy;
use File::Spec;
use Getopt::Long;
use ReseqTrack::Tools::QAandCntFastq;
use ReseqTrack::Tools::Counts;


$| = 1; 

my $start_run = time();

my (
	$run_id,
	$output_dir,
#	$log_dir,
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $help,
    $type,
    $run_dir,
    $len_limit,
   );

my  $rsync_options = '-avc --rsh=/usr/bin/rsh';
my  $load = 0;
my  $no_update = 0;
my  $backfill_stats = 0;
my  $host_name = "1000genomes.ebi.ac.uk";
my  $remote = 0;
my  $new_type = "FILTERED_FASTQ";
my $clobber = 0;

&GetOptions(
  'run_id=s'		=>	\$run_id,
  'type=s'			=>	\$type,
  'output_dir:s'	=>	\$output_dir,
#  'log_dir=s' 		=>	\$log_dir,
  'run_dir:s' 		=> 	\$run_dir,
  'rsync_options=s' =>	\$rsync_options,
  'dbhost=s'     => \$dbhost,
  'dbname=s'     => \$dbname,
  'dbuser=s'     => \$dbuser,
  'dbpass=s'     => \$dbpass,
  'dbport=s'     => \$dbport,
  'help!'		 => \$help,
  'load!' 	 	 => \$load,
  'host_name=s'	 => \$host_name,
  'remote!'     => \$remote,
  'backfill_stats!' => \$backfill_stats,
  'no_update!'			=> \$no_update,
  'len_limit:i'		=> \$len_limit,
  'clobber!' => \$clobber,
);
if ($help) {
	help_info();
}

if (!$run_id || !$type ) {
	throw( "please check if run_id and type are entered at the command line\n");
}
#if ( !$log_dir ) {   
#	throw("please check if log_dir is entered at the command line.\n");
#}	
	
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

# query information from Collection table    
my $ca = $db->get_CollectionAdaptor;
my $collection = $ca->fetch_by_name_and_type($run_id, $type);
throw("No collection object is found; please check the run_id and type\n") if (!$collection);
my $others = $collection->others; #return a reference of an array of FileAdaptor objects. a FileAdaptor object contains all info in a row in the File table
throw("No File object is found for run_id $run_id\n") if (!$others);
my $new_collection = $ca->fetch_by_name_and_type($run_id, $new_type);
if($new_collection){
  print $run_id." already seems to have a ".$new_type." collection with \n";
  foreach my $other(@{$new_collection->others}){
    print $other->name."\n";
  }
  exit(0) unless($clobber);
}
# Query info from meta_run_info table
my $rmi_a = $db->get_RunMetaInfoAdaptor;
my $rmi = $rmi_a->fetch_by_run_id($run_id);
if($rmi->status eq 'suppressed'){
  print $run_id." is now suppressed in the archive, skipping the run\n";
  exit(0);
}
#check if a single ended run has multiple files
my ($mate1, $mate2, $frag) = assign_file_objects($others);
if($rmi->library_layout eq 'SINGLE'){
  if(!$frag || ($mate1 || $mate2)){
    print STDERR "There is a problem for ".$rmi->run_id." is has the wrong files\n";
    print STDERR "There are ".@$others." files\n";
    foreach my $file(@$others){
      print STDERR $file->name."\n";
    }
    print STDERR $rmi->run_id." has library layout ".$rmi->library_layout."\n";
    throw("ReseqTrack/scripts/qc/run_fastq_qc.pl can't run filtered on ".
          $collection->name." it has the wrong number of files");
  }
}else{
  unless($mate1 && $mate2){
    print STDERR "There is a problem for ".$rmi->run_id." is has the wrong files\n";
    foreach my $file(@$others){
      print STDERR $file->name."\n";
    }
    throw("ReseqTrack/scripts/qc/run_fastq_qc.pl can't run filtered on ".
          $collection->name." it has the wrong number of files"); 
  }
}
#print "Setting disconnect when inactive to 1\n";
$db->dbc->disconnect_when_inactive(1);

# set up output directory if -output_dir option is entered
my $full_output_dir;
if($output_dir){
  my $individual = $rmi->sample_name;
  $full_output_dir = $output_dir."/".$individual."/sequence_read";
  $full_output_dir =~ s/\/\//\//;
}

#my $log_file_path = $log_dir . "/" . $run_id . ".log";
	
my %files; 	#### key is the file path after it is copied to a run dir if -run_dir option is specified; 
			#### when -run_dir is not specified, key is the original file path.
			#### value is the original file object
foreach my $other (@$others) {
  my $filePath = $other->name;
  #print "input file name: $filePath\n";
  
  #### if output_dir is not given on command line, write filtered fastq files one level above the input file, in a dir called "filtered_sequence"
  if (!$full_output_dir) {
    my ($base_filename, $dir, $suffix) = fileparse($filePath);	
    my @array = split(/\//, $dir);
    pop @array;
    my $upper_dir = File::Spec->catdir(@array);	
    #print "upper dir is $upper_dir\n";	
    $full_output_dir = $upper_dir . "/sequence_read";
  }	
  mkpath($full_output_dir) unless (-d $full_output_dir);
  unless(-d $full_output_dir){
    throw("Can't run without ".$full_output_dir." existing");
  }
  if ($run_dir) {
    my $movedFastqName = copy_fastq($filePath);		
    $files{$movedFastqName} = $other;
  }
  else {
    $files{$filePath} = $other;		
  }	
}

#### $outfiles is a referench to a hash
#### key is filtered fastq file pathes, 
#### value is an array, 1st element is Counts object for the file (stats), second element is the original file object
#my $outfiles = qa_and_cnt(\%files, $rmi, \$full_output_dir, \$log_file_path, $len_limit);
my $outfiles = qa_and_cnt(\%files, $rmi->instrument_platform, \$full_output_dir, $len_limit, $run_id);

throw("ERROR: no qa output files exist\n") if ( !%$outfiles);

my (	%outfile_md5,  #### These three hashes, key is the zipped filtered fastq file paths; value is md5, stats obj and original fiastq file obj, sperately
		%outfile_stats, 
		%outfile_ori_obj,
	);
foreach my $outfile (keys %$outfiles) {
	my $size = -s $outfile;
	if ($size == 0) { # newly added
		#print "size 0 files $outfile\n";
		`rm $outfile`;
		next;
	}	
	`gzip -f $outfile`;
	my $zipped_outfile = $outfile . ".gz";
	my $zipped_outfile_md5 = run_md5($zipped_outfile); #generare md5 for zipped files
	
	print "file path: $outfile\t$zipped_outfile_md5\tsize $size\n";

	$outfile_md5{$zipped_outfile}=$zipped_outfile_md5;
	
	my $tmp_array = $$outfiles{$outfile}; 				#anonymous array
	$outfile_stats{$zipped_outfile}=$$tmp_array[0]; 	#hash, key is file path, value is the Counts object for the file
	$outfile_ori_obj{$zipped_outfile}=$$tmp_array[1]; 	#hash, key is file path, value is the oringinal unfiltered fastq file object
}	

load_files(\%outfile_md5, $new_type, \%outfile_stats, \%outfile_ori_obj) if ($load);

if ($run_dir) {	
	foreach my $f (keys %files) {
		`rm -f $f`;
	}
}		

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";

########################## SUBS #################################################
sub load_files {
	
	#### This function load file objects into tables file, collection_group and history  												####
	#### It stores base and read counts into the statistics table for the filtered fastq files 											####
	#### It stores base and read countd of the ORIGINAL fastq files in the statistics table when -backfill_stats is used 				####
	#### If -update option is entered and the script is re-run for the same input file, the filtered fastq file will be overwritten and ####
	#### history table will have an update record, statistics table will be updated as well 											####
	 
	my ($filePath_md5, $new_file_type, $stats, $ori_file_obj) = @_;

	$db->dbc->disconnect_when_inactive(0);
	
	my $ha = $db->get_HostAdaptor;
	my $host = $ha->fetch_by_name($host_name);
	if(!$host){
	  $host = ReseqTrack::Host->new
	      (
	       -name => $host_name,
	       -remote => $remote
	      );
	}
	
	my @objects;
	my $fa = $db->get_FileAdaptor;
	
	foreach my $path(keys(%$filePath_md5)){
	  #print "Inside load: path is: $path\n";
	  my $files = create_objects_from_path_list([$path], $new_file_type, $host); #$files is a reference to an array of file objects (in this case, only one element)
	  my $file = $files->[0]; #to get the first and only file object
	  my $md5 = $filePath_md5->{$file->name};
	  if(!$md5){
	    warning("Don't have a md5 for ".$file->filename);
	  }
	  $file->md5($md5); #set md5

	  #get read and base counts from the Counts object	  
	  my $read_count = $$stats{$path}->filt_read_cnt;
	  my $base_count = $$stats{$path}->filt_base_cnt;
	  my $read_count_b4qc = $$stats{$path}->read_cnt;
	  my $base_count_b4qc = $$stats{$path}->base_cnt;
	 
	  ### skip loading files if the filtered files are empty, still fill in the stats forthe original files though
	  if (!$read_count && !$base_count) {
	      warning("the filtered file $path is empty\n");
	      goto BACKFILL_STATS;
	  } 
	     
	  print "read cnt before QA: " . $read_count_b4qc . "\n" if ($read_count_b4qc);
	  print "base cnt before QA: " . $base_count_b4qc . "\n" if ($base_count_b4qc);
	  print "read cnt after  QA: " . $read_count . "\n" if ($read_count);
	  print "base cnt after  QA: " . $base_count . "\n" if ($base_count);
	  
	  #create statistics for file object
      my $read_stats = create_statistic_for_object($file, 'read_count', $read_count);
      my $base_stats = create_statistic_for_object($file, 'base_count', $base_count);
      
      #give file object statistics objects
      $file->statistics($read_stats);
      $file->statistics($base_stats);
     
     #### Everytime when the QA script is run with a fresh input file, a history object would created for the output filtered fastq file and inserted in the history table
	 #### When the filtered fastq file object already exist in the file table (identical basename of files), 
	 #### update the file object, use the old dbID in file table
	 #### get the exisitng file's history object, change the comment field and store it as a new row in the history table
	 
	 my $history;
	 my $hist_a = $db->get_HistoryAdaptor();
	 my $possible_stored_file = $fa->fetch_by_filename($file->filename); #filename function gets the basename of the file; name function gets the whole path 	
	 
	 if($possible_stored_file && @$possible_stored_file >0 ){
	     if ( $no_update ) {
	     	throw("Filtered fastq file already exist somewhere in the system, cannot continue if the -no_update option is set");
	     } 
	     if(@$possible_stored_file == 1) {
	        #### get the possible stored object
	        my $stored_file = $possible_stored_file->[0];	
	        #print "matching filename stored in db: " . $stored_file->name . "\n";	
	        
	        #### compare new file and the exisitng file to see what the difference between the files
	        my $comment = calculate_comment($file, $stored_file); 
			
			$history = ReseqTrack::History->new(
				-other_id => $stored_file->dbID,
				-table_name => 'file',
				-comment => $comment, 
	 		);
		
			#### update the history object and assign it to the new file
	 		$file->history($history);
			
			#### modify the file object and store it in the file table with the update function
			$file->dbID($stored_file->dbID);
	 		$file->created($stored_file->created);
	 		$fa->update($file,1); #1 is withdraw indicator for overwrite the old file in file table 
	     }
	     elsif ( @$possible_stored_file > 1) {
	     	 warning("There are multiple files in the database with the filename ".
	                  $file->filename." Don't know what to do now");
		 }
	 }
	 #### If the file object does not exist in db, store it and create a history object for it
	 else{
		$fa->store($file);	
		$history = ReseqTrack::History->new(
			-other_id => $file->dbID,
			-table_name => 'file',
			-comment => "Created new filtered fastq file", 
	 	);
	 	$hist_a->store($history);
	 	$file->history($history);
	 	#print "now file has a history: ". $file->history->[0]->comment . "\n";	
	}
	   
	 push(@objects, $file);
	
	 BACKFILL_STATS:	 
	
	 #### when $backfill_stats is used, populate statistics table with counts for the original, unfiltered fastq files	
	 if ($backfill_stats) {	  	  
		#### get the file object for the original unfiltered file
		my $unfilt_file_obj = $$ori_file_obj{$path};
		
		#print "\nI am inside back fill stats block\n";
		if ( $base_count_b4qc && $read_count_b4qc && $unfilt_file_obj ) { #files like ERR000587.filt.fastq.gz won't always have an original unfiltered 
																			#fastq file associated with it	 
			my $unfilt_file_read_stats = create_statistic_for_object($unfilt_file_obj, 'read_count', $read_count_b4qc);
			my $unfilt_file_base_stats = create_statistic_for_object($unfilt_file_obj, 'base_count', $base_count_b4qc);	  
			
			$unfilt_file_obj->statistics($unfilt_file_read_stats);
			$unfilt_file_obj->statistics($unfilt_file_base_stats);

			if (!$unfilt_file_obj->statistics->[0]->dbID || !$unfilt_file_obj->statistics->[1]->dbID) { # if the stas object does not have a dbID
			    $fa->store_statistics($unfilt_file_obj);
			}
			else {
			    $fa->update_statistics($unfilt_file_obj); 
			    #print "unfilt_file_obj stats obj: " . $unfilt_file_obj->statistics->[0]->attribute_value . "\n"; 
			  	#print "unfilt_file_obj stats obj: " . $unfilt_file_obj->statistics->[1]->attribute_value . "\n";
			  	#print "unfilt_file_obj stats obj dbID: " . $unfilt_file_obj->statistics->[0]->dbID . "\n";
			  	#print "unfilt_file_obj stats obj dbID: " . $unfilt_file_obj->statistics->[1]->dbID . "\n";
			  	#print "unfilt_file_obj stats obj other id: " . $unfilt_file_obj->statistics->[0]->other_id . "\n";
			  	#print "unfilt_file_obj stats obj other id: " . $unfilt_file_obj->statistics->[1]->other_id . "\n";
			}          
		}
		#print "I just come out of the backfill block\n\n";     
	  }
	} # the end of foreach path loop
	if(!@objects || @objects == 0){
          throw("Failed to produce any fastq files of ".$new_file_type." for ".$run_id);
        }
	my $collection =  ReseqTrack::Collection->new(
	  -name => $run_id,
	  -others => \@objects,
	  -type => $new_file_type,
	  -table_name => 'file',
	);
	
	my $ca = $db->get_CollectionAdaptor; 
	$ca->store($collection) if($load);
	
	return 1;
}	

sub copy_fastq {
	my ($fPath) = @_;
	my $fName = basename( $fPath);
	my $rsync = "rsync ".$rsync_options." ".$fPath." ".$run_dir."/.";
 	eval{
    	system($rsync);
  	};
  	if($@){
    	throw("Failed to run rsync ".$rsync." $@");
  	}
  	my $fileName = $run_dir . "/" . $fName;

 	return ($fileName);
		
}	

sub help_info {

 exec('perldoc', $0);

}

#############################################################################################
=pod

=head1 NAME

 perl -w /homes/zheng/reseq-personal/zheng/bin/runQC.pl

 Required arguments:
   
	-dbhost, the name of the mysql-host
	-dbname, the name of the mysql database
	-dbuser, the name of the mysql user
	-dbpass, the database password if appropriate
	-dbport, the port the mysql instance is running on, this defaults to 4197 the standard
             port for mysql-g1kdcc.ebi.ac.uk

	-run_id,
	-type,

 Arguments with default settings:
   
	-load = 0			use -load when you want to load the end results into the tracking db 	
	-backfill_stats = 0	use -backfill_stats when you want to store base and read counts for original unfiltered fastq files in the statistics table
	
	-host_name = "1000genomes.ebi.ac.uk"
	-remote = 0
	-rsync_options = '-avc --rsh=/usr/bin/rsh'

 OPtional arguments:
 
    -help			This makes the script print out its options
    -run_dir		to specify the run directory, fastq files will be copied there and all QA will take place there
    -len_limit		as a default, for Illumin and Solid, QA checks the first 25bp; for 454, it is 30bp. This option allows you to change that cutoff length 
	-output_dir		specify where the filtered fastq file should be output; default is in a dir called "filtered_sequence" above the input dir
	-no_update		No update will be performed if -no_update option is set.  If the flag is NOT set, if the same input file is re-run through the QA pipeline, 
					the output filtered fastq files will be overwritten and statistics will be overwritten; a new row will be inserted in the history table.   
	    	
=head1 SYNOPSIS

This script take a run_id, find all fastq files associated with it from the tracking db, then QA the fastq files. Below is a list of things
the QA process checks:

Syntax Checks:

	-Each header line begins with @
	-The third line always starts with a +
	-There are four lines in each entry (implied by the above two rules)
	-On line3, if there is a name following the + sign, the name has to match the one found in line1
	-The sequence and quality lines are the same length
	-When there are paired end files that both the _1 and _2 files have the same number of reads in them. 
	-For SOLID colourspace fastq, each read starts with a base before going into the numbers which represent colourspace

Sequence Checks:

	-Read is longer than 35bp for Solexa, 25bp for Solid, and 30 bp for 454
	-Read does not contain any N's in the first 25, 30 or 35bp
	-Quality values are all 2 or higher in the first 25bp, 30bp or 35bp
	-The reads contain more than one type of base in the first 25, 30, or 35bp
	-Read does not contain more than 50% Ns in its whole length

The output of the script are the followings:

	- fatsq files of filtered sequences; if one of the paired reads failed, the good one is moved to the fragment.fastq file 
	- the location of the above files are loaded in the tracking db, if -load is used
	- the base count and read count of the fatsq file AFTER QA is stored in the statistics table if -load is used;
	- the base count and read count of the fatsq file BEFORE QA is stored in the statistics table if -backfill_stats is used
	
=head1 Example:

 perl -w ~/ReseqTrack/scripts/qc/run_fastq_qc.pl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -dbname zheng_fastq_test -run_id ERR000587 -load -type FASTQ -output_dir /tmp > out

 perl -w ~/reseq-personal/zheng/bin/runQC.pl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -dbname zheng_fastq_test 
 -run_id ERR000587 -load -backfill_stats -type ARCHIVE_FASTQ -output_dir /tmp  > out

	ERR000587 is a test case for ILLUMINA with two files
	ERR000591 is a test case for ILLUMINA with three files
	ERR001623 is a test case for Solid with one file
	SRR014948 is a test case for 454 with two files
 
=head1 Author:

 Holly Zheng Bradley Oct. 2009

=cut
