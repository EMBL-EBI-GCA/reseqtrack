
#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FilterFastq;
use ReseqTrack::Tools::GeneralUtils;
use File::Path;
use ReseqTrack::Tools::Loader::File;
use ReseqTrack::Tools::Loader::Archive;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Tools::StatisticsUtils;
use ReseqTrack::Collection;
use Getopt::Long;

$| = 1;

my ($dbhost, $dbuser, $dbpass, $dbname, $dbport);
my $run_id;
my $collection_type;
my $new_collection_type = 'FILTERED_FASTQ';
my $clobber;
my $final_dir_name = 'sequence_read';
my $output_dir;
my $host_name = '1000genomes.ebi.ac.uk';
my $min_length = 25;
my $min_qual = 2;
my $max_percent_n = 0.5;
my $use_platform_default_length = 1;
my $help;

&GetOptions(
	    'dbhost=s'     => \$dbhost,
	    'dbname=s'     => \$dbname,
	    'dbuser=s'     => \$dbuser,
	    'dbpass=s'     => \$dbpass,
	    'dbport=s'     => \$dbport,
	    'run_id=s' => \$run_id,
	    'output_dir=s' => \$output_dir,
	    'host_name=s' => => \$host_name,
	    'collection_type=s' => \$collection_type,
	    'new_collection_type=s' => \$new_collection_type,
	    'final_dir_name=s' => \$final_dir_name,
	    'clobber!' => \$clobber,
	    'min_length=s' => \$min_length,
	    'min_qual=s' => \$min_qual,
	    'max_percent_n=s' => \$max_percent_n,
	    'use_platform_default_length!' => \$use_platform_default_length,
	    'help!' => \$help,
	   );

if($help){
  useage();
}

#Creating database connection
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

unless($db){
  throw("Can't run without a database");
}

#Need run id and collection type to get required files and meta info
unless($run_id && $collection_type){
  throw("Can't run without a run id and a collection type");
}

my $ca = $db->get_CollectionAdaptor;
my $collection = $ca->fetch_by_name_and_type($run_id, $collection_type);

throw("Failed to find a collection for ".$run_id." ".$collection_type) 
  unless($collection);
#Check if output collection already exists
my $new_collection = $ca->fetch_by_name_and_type($run_id, $new_collection_type);
if($new_collection){
  print $new_collection->name." ".$new_collection_type." already exists\n";
  unless($clobber){
    print "Exiting process\n";
    exit(0);
  }
}
#get input info
my $others = $collection->others;
my $rmia = $db->get_RunMetaInfoAdaptor;
my $rmi = $rmia->fetch_by_run_id($run_id);
throw("Can't run without meta info for ".$run_id) unless($rmi);
if($use_platform_default_length){
  $min_length = 35 if($rmi->instrument_platform eq 'ILLUMINA');
  $min_length = 30 if($rmi->instrument_platform eq 'LS454');
}
#Skip non public meta info
if($rmi->status ne 'public'){
  print $run_id." is now ".$rmi->status." in the archive so skipping\n";
  exit(0);
}
#Check you have the correct input files
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


$db->dbc->disconnect_when_inactive(1);
#Create output directory
my $full_output_dir = $output_dir."/".$rmi->sample_name."/".$final_dir_name;
mkpath($full_output_dir) unless(-d $full_output_dir);
throw("Full output dir ".$full_output_dir." does not exist") unless(-d $full_output_dir);
#Define if the sequence is colourspace
my $is_colorspace = 0;
$is_colorspace = 1 if($rmi->instrument_platform eq 'ABI_SOLID');
#print "Have ".@$others." from ".$collection->name."\n";

my $mate1_name = $mate1->name if($mate1);
my $mate2_name = $mate2->name if($mate2);
my $frag_name = $frag->name if($frag);
#create FilterFastq object
my $filter_fastq = ReseqTrack::Tools::FilterFastq->new
  (
   -mate1 => $mate1_name,
   -mate2 => $mate2_name,
   -frag => $frag_name,
   -run_id => $rmi->run_id,
   -output_dir => $full_output_dir,
   -is_colorspace => $is_colorspace,
   -min_length => $min_length,
   -min_qual => $min_qual,
   -max_percent_n => $max_percent_n,
   -clobber => $clobber,
  );
#Run Filtering
my $files = $filter_fastq->filter_files;
#Check the files returned all exist
foreach my $file(@$files){
  unless($file =~ /\.fastq\.gz$/){
    throw("Not sure what to do with ".$file." doesn't match expect file extension");
  }
}
exit(0) unless($files && @$files >= 1);
#create a File loader object for the files
my $loader = ReseqTrack::Tools::Loader::File->new
  (
   -file => $files,
   -do_md5 => 1,
   -hostname => $host_name,
   -db => $db,
   -assign_types => 0,
   -check_types => 0,
   -type => $new_collection_type,
   -update_existing => 1,
  );
#Load the files into the database
$loader->process_input();
$loader->create_objects();
$loader->sanity_check_objects();
my $file_objects = $loader->load_objects();
#create filtered collection
if(!$file_objects || @$file_objects == 0){
  throw("Have no file objects to create collection with");
}
#Create a collection
my $filtered_collection = ReseqTrack::Collection->new(
						      -name => $rmi->run_id,
						      -others => $file_objects,
						      -type => $new_collection_type,
						      -table_name => 'file',
						     );
#Store the collection
$db->get_CollectionAdaptor->store($filtered_collection);
my ($filt_m1, $filt_m2, $filt_f) = assign_file_objects($file_objects);
#Update statistics
my @objects_to_update;
if($mate1 && $mate2){
  my $mate1_rc = create_statistic_for_object
    ($mate1, 'read_count', $filter_fastq->unfiltered_mate1_readcount);
  $mate1->statistics($mate1_rc);
  my $mate1_bc = create_statistic_for_object
    ($mate1, 'base_count', $filter_fastq->unfiltered_mate1_basecount);
  $mate1->statistics($mate1_rc);
  my $mate2_rc = create_statistic_for_object
    ($mate2, 'read_count', $filter_fastq->unfiltered_mate2_readcount);
  $mate2->statistics($mate2_rc);
  my $mate2_bc = create_statistic_for_object
    ($mate2, 'base_count', $filter_fastq->unfiltered_mate2_basecount);
  $mate2->statistics($mate2_rc);
  push(@objects_to_update, $mate1, $mate2);
}else{
  if(($mate1 && !$mate2) || (!$mate1 && $mate2)){
    throw("There seems to be an issue we only have either mate 1".$mate1." or ".
	  " mate2 ".$mate2." not both");
  }
}
if($frag){
  my $frag_rc = create_statistic_for_object
    ($frag, 'read_count', $filter_fastq->unfiltered_frag_readcount);
  $frag->statistics($frag_rc);
  my $frag_bc = create_statistic_for_object
    ($frag, 'base_count', $filter_fastq->unfiltered_frag_basecount);
  $frag->statistics($frag_bc);
  push(@objects_to_update, $frag);
}
if($filt_m1 && $filt_m2){
  print "Creating statistic for ".$filt_m1->name."\n";
  my $filt_m1_rc = create_statistic_for_object
    ($filt_m1, 'read_count', $filter_fastq->filtered_mate1_readcount);
  $filt_m1->statistics($filt_m1_rc);
  my $filt_m1_bc = create_statistic_for_object
    ($filt_m1, 'base_count', $filter_fastq->filtered_mate1_basecount);
  $filt_m1->statistics($filt_m1_bc);
  my $filt_m2_rc = create_statistic_for_object
    ($filt_m2, 'read_count', $filter_fastq->filtered_mate2_readcount);
  $filt_m2->statistics($filt_m2_rc);
  my $filt_m2_bc = create_statistic_for_object
    ($filt_m2, 'base_count', $filter_fastq->filtered_mate2_basecount);
  $filt_m2->statistics($filt_m2_bc);
  push(@objects_to_update, $filt_m1, $filt_m2);
}else{
  if(($filt_m1 && !$filt_m2) || (!$filt_m1 && $filt_m2)){
    throw("There seems to be an issue we only have either filt_m 1".$filt_m1." or ".
	  " filt_m2 ".$filt_m2." not both");
  }
}
<<<<<<< .mine

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
	
	print "Path is $zipped_outfile\n";
	print "read_count = " . $$tmp_array[0]->filt_read_cnt ."\n";
        print "base_count = " . $$tmp_array[0]->filt_base_cnt . "\n";
	print "unfiltered_read_count = " . $$tmp_array[0]->read_cnt ."\n";
        print "unfiltered_base_count = " . $$tmp_array[0]->base_cnt . "\n";
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
          warning("Failed to produce any fastq files of ".$new_file_type." for ".$run_id);
	  return;
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

=======
if($filt_f){
  my $frag_rc = create_statistic_for_object
    ($filt_f, 'read_count', $filter_fastq->filtered_frag_readcount);
  $filt_f->statistics($frag_rc);
  my $frag_bc = create_statistic_for_object
    ($filt_f, 'base_count', $filter_fastq->filtered_frag_basecount);
  $filt_f->statistics($frag_bc);
  push(@objects_to_update, $filt_f);
>>>>>>> .r270
}
#Store all the statistic objects
my $fa = $db->get_FileAdaptor;
foreach my $file(@objects_to_update){
  $fa->store_statistics($file, 1);
}

my $archiver = ReseqTrack::Tools::Loader::Archive->new
  (
   -file      => $files,
   -action    => 'archive',
   -priority  => 90,
   -no_lock   => 1,
   -db => $db,
  );
$archiver->process_input();
$archiver->sanity_check_objects();
$archiver->archive_objects();

