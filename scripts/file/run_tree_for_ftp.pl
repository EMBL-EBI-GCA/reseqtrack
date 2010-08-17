#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor; 
use ReseqTrack::File;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::ArchiveUtils;
use File::Basename;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $run = 0;
my $dir_to_tree = '/nfs/1000g-archive/vol1/ftp';
my $output_path = '/nfs/1000g-work/G1K/archive_staging/ftp/current.tree';


&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'dir_to_tree=s' => \$dir_to_tree,
  'output_path=s' => \$output_path,
    );



my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

my $fa = $db->get_FileAdaptor;
#fetch old tree object
my $files = $fa->fetch_by_filename("current.tree");
if(!$files || @$files == 0){
  throw("Can't find the current.tree file in ".$dbname);
}elsif(@$files >= 2){
  throw("Seem to have ".@$files." versions of the current.tree file don't know ".
        "what to do");
}

my $old_tree = $files->[0];
throw("Old tree ".$old_tree->name." doesn't have a dbID in ".$dbname) unless($old_tree->dbID);
#run tree
dump_dirtree_summary($dir_to_tree, $output_path, undef, $fa);
#create new file object
my ($filename, $location_path) = fileparse($output_path);
my $new_tree =  create_object_from_path($output_path, $old_tree->type, $old_tree->host);
my $md5 = run_md5($new_tree->full_path);
$new_tree->md5($md5);

my $archive_action_adaptor = $db->get_ArchiveActionAdaptor;
my $archive_location_adaptor = $db->get_ArchiveLocationAdaptor;
my $aa = $db->get_ArchiveAdaptor;

my $archives = $aa->fetch_all;
cleanup_archive($archives, $db, 0);

my $action = $archive_action_adaptor->fetch_by_action('replace');
my $archive_location = $archive_location_adaptor->fetch_by_archive_location_name
    ('staging');

$aa->delete_archive_lock;


if($new_tree->md5 ne $old_tree->md5){
  #store new object
  throw("Don't have dbID for ".$old_tree->name) unless($old_tree->dbID);
  my $history = create_history($new_tree, $old_tree);
  $new_tree->history($history);
  $fa->store($new_tree, 1);
  #create archive
  my $archive_root = $archive_location->location;
  throw("Can't archive ".$new_tree->full_path." it doesn't live in ".$archive_location->location) unless($new_tree->full_path =~ /$archive_root/);
  my $archive = create_archive_from_objects($new_tree, $action, $archive_location);
  $archive->priority(90);
  $aa->store($archive);
  my $true = 1;
  while($true){
    my $archives = $aa->fetch_all;
    my @archives;
    foreach my $archive(@$archives){
      push(@archives, $archive) unless($archive->fire_exit_code);
    }
    $true = 0 if(@archives == 0);
    cleanup_archive($archives, $db, 0);
    sleep(10);
  }
 $aa->delete_archive_lock;	
}else{
  print STDERR "The tree files are the same\n";
  unlink $new_tree->full_path;
}


