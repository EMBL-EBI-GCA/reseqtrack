#!/sw/arch/bin/perl -w

use strict;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::ArchiveUtils;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::Exception;
use File::Basename;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $move_list;
my $clobber = 0;
my $run = 0;
my $priority = 90;
my $from;
my $to;

&GetOptions( 
	    'dbhost=s'      => \$dbhost,
	    'dbname=s'      => \$dbname,
	    'dbuser=s'      => \$dbuser,
	    'dbpass=s'      => \$dbpass,
	    'dbport=s'      => \$dbport,
	    'move_list=s' => \$move_list,
	    'from=s' => \$from,
	    'to=s' => \$to,
	    'priority=s' => \$priority,
	    'clobber!' => \$clobber,
  'run!' => \$run,
    );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

my $archive_action_adaptor = $db->get_ArchiveActionAdaptor;
my $archive_location_adaptor = $db->get_ArchiveLocationAdaptor;
my $aa =  $db->get_ArchiveAdaptor;

my $action = $archive_action_adaptor->fetch_by_action('move_within_volume');
throw("Failed to get action for move_within_volume") unless($action);
my $location = $archive_location_adaptor->fetch_by_archive_location_name('archive');

my $hash;
if($move_list){
  $hash = parse_movelist($move_list);
}elsif($from && $to){
  $hash->{$from} = $to;
}else{
  throw("Need to give move_archive_files.pl either a -move_list or a -from and ".
	"a -to value");
}
my $fa = $db->get_FileAdaptor;

my @file_objects;
my $location_root = $location->location;
foreach my $path(keys(%$hash)){
  my $new_path = $hash->{$path};
  throw($path." doesn't exist can't move it") unless(-e $path);
  if(!$clobber && -e $new_path){
    warning($new_path." exists can't move ".$path." on top of it");
    next;
  }
  throw("Either ".$path." or ".$new_path." isn't in the archive can't move ".
        "them ".$location->location) unless($path =~ /$location_root/ && 
                        $new_path =~ /$location_root/);
  my $name = basename($path);
  my $object = $fa->fetch_by_name($path);
  throw("Failed to retrive a file object for ".$path) unless($object);
  push(@file_objects, $object);
}

print "Have ".@file_objects." to move \n";
my @archives;
foreach my $file(@file_objects){
  my $new_path = $hash->{$file->full_path};
  throw ("Don't seem to have a new path for ".$file->full_path) unless($new_path);
  throw("File doesn't exist ".$file->name) unless(-e $file->name);
  my $archive = create_archive_from_objects($file, $action, $location, $new_path);
  push(@archives, $archive);
}

print "Have ".@archives." archive objects\n";

foreach my $archive(@archives){
  $archive->priority($priority);
  $aa->store($archive) if($run);
}


$aa->delete_archive_lock;
