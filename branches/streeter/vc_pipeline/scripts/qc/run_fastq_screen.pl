use strict;
use warnings;

use ReseqTrack::Tools::SequenceIndexUtils qw (assign_files);
use ReseqTrack::Tools::StatisticsUtils;
use ReseqTrack::Tools::QC::FastQScreen;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::StatisticsUtils;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Tools::Loader::File;
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use Getopt::Long;
use File::Copy;
use File::Path qw(make_path);

use Data::Dumper;

my ($dbhost, $dbuser, $dbpass, $dbname, $dbport);
my $collection_name;
my $collection_type;
my $output_dir;
my $host_name = '1000genomes.ebi.ac.uk';
my $fastqscreen_path = 'fastqscreen';
my $conf_file;
my $directory_layout;

&GetOptions(
	'dbhost=s'     => \$dbhost,
	'dbname=s'     => \$dbname,
	'dbuser=s'     => \$dbuser,
	'dbpass=s'     => \$dbpass,
	'dbport=s'     => \$dbport,
	'collection_name=s' => \$collection_name,
	'output_dir=s' => \$output_dir,
	'host_name=s' => => \$host_name,
	'collection_type=s' => \$collection_type,
	'fastqscreen_path=s' => \$fastqscreen_path, 
	'conf_file=s' => \$conf_file,
	'directory_layout=s' => \$directory_layout,
);
	

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
unless($collection_name && $collection_type){
  throw("Can't run without a collection name and a collection type");
}

my $ca = $db->get_CollectionAdaptor;
my $file_adaptor = $db->get_FileAdaptor;
my $collection = $ca->fetch_by_name_and_type($collection_name, $collection_type);

throw("Failed to find a collection for $collection_name $collection_type from ".$ca->dbc->dbname) 
  unless($collection);

throw("output_dir must be specified") unless ($output_dir);

my $rmia = $db->get_RunMetaInfoAdaptor;
my $run_meta_info = $rmia->fetch_by_run_id($collection_name);

if ($run_meta_info && $directory_layout) {
	$output_dir = create_directory_path($run_meta_info, $directory_layout, $output_dir);
}


if (! -d $output_dir) {
	make_path($output_dir) or throw("Could not create output dir $output_dir");
}


my @files = assign_files($collection->others);

if (defined $files[0] && defined $files[1]) {
	@files = ($files[0],$files[1]);
}
else {
	@files = ($files[2]);
}



my @input_files = map {$_->name} @files;


my $fastqscreen = ReseqTrack::Tools::QC::FastQScreen->new(
	-program => $fastqscreen_path,
	-keep_text => 1,
	-keep_graph => 1,
	-working_dir => $output_dir,
	-input_files => \@input_files,
	-conf_file => $conf_file,
);

$fastqscreen->run;


my $graph_destination = "$output_dir/${collection_name}_fastqscreen.png";
my $report_destination =  "$output_dir/${collection_name}_fastqscreen.txt";
move($fastqscreen->graph_path,$graph_destination);
move($fastqscreen->text_path,$report_destination);


#create a File loader object for the files
my $loader = ReseqTrack::Tools::Loader::File->new
  (
   -file => [$graph_destination,$report_destination],
   -do_md5 => 1,
   -hostname => $host_name,
   -db => $db,
   -assign_types => 0,
   -check_types => 0,
   -type => 'FASTQSCREEN',
   -update_existing => 1,
  );
#Load the files into the database
$loader->process_input();
$loader->create_objects();
$loader->sanity_check_objects();
my $file_objects = $loader->load_objects();

my $output_collection = ReseqTrack::Collection->new(
  -name => $collection_name,
  -others => $file_objects,
  -type => 'FASTQSCREEN',
  -table_name => 'file',
);

$db->get_CollectionAdaptor->store($output_collection);

