
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

throw("Failed to find a collection for ".$run_id." ".$collection_type." from ".$ca->dbc->dbname) 
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
my ($mate1, $mate2, $frag) = assign_files($others);
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
#$db->dbc->disconnect_when_inactive(0);
foreach my $file(@$files){
  unless($file =~ /\.fastq\.gz$/){
    throw("Not sure what to do with ".$file." doesn't match expect file extension");
  }
}
exit(0) unless($files && @$files >= 1);
$db->dbc->disconnect_when_inactive(0);
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
#$db->get_CollectionAdaptor->store($filtered_collection);
my ($filt_m1, $filt_m2, $filt_f) = assign_files($file_objects);
#Update statistics
my @objects_to_update;
if($mate1 && $mate2){
  my $mate1_rc = create_statistic_for_object
    ($mate1, 'read_count', $filter_fastq->unfiltered_mate1_readcount);
  $mate1->statistics($mate1_rc);
  my $mate1_bc = create_statistic_for_object
    ($mate1, 'base_count', $filter_fastq->unfiltered_mate1_basecount);
  $mate1->statistics($mate1_bc);
  my $mate2_rc = create_statistic_for_object
    ($mate2, 'read_count', $filter_fastq->unfiltered_mate2_readcount);
  $mate2->statistics($mate2_rc);
  my $mate2_bc = create_statistic_for_object
    ($mate2, 'base_count', $filter_fastq->unfiltered_mate2_basecount);
  $mate2->statistics($mate2_bc);
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
if($filt_f){
  my $frag_rc = create_statistic_for_object
    ($filt_f, 'read_count', $filter_fastq->filtered_frag_readcount);
  $filt_f->statistics($frag_rc);
  my $frag_bc = create_statistic_for_object
    ($filt_f, 'base_count', $filter_fastq->filtered_frag_basecount);
  $filt_f->statistics($frag_bc);
  push(@objects_to_update, $filt_f);
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

