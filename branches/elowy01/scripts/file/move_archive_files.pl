#!/usr/bin/env perl

use strict;
use warnings;

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
my $help;

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
	    'help!' => \$help,
    );

if($help){
  useage();
}

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

###print "Have ".@file_objects." to move \n";
my @archives;
foreach my $file(@file_objects){
  my $new_path = $hash->{$file->full_path};
  throw ("Don't seem to have a new path for ".$file->full_path) unless($new_path);
  throw("File doesn't exist ".$file->name) unless(-e $file->name);
  my $archive = create_archive_from_objects($file, $action, $location, $new_path);
  #print "moving *".$file->name."* to *".$new_path."*\n";
  push(@archives, $archive);
}

print "Have ".@archives." archive objects\n";

foreach my $archive(@archives){
  $archive->priority($priority);
  $aa->store($archive) if($run);
}


$aa->delete_archive_lock;


print "\nMessage:\n";
print " move_archive_files.pl does not automatically\n";
print " run 'cleanup_archive.pl' on archive table.\n";

=pod

=head1 NAME

ReseqTrack/scripts/files/move_archive_files.pl

=head1 SYNOPSIS

This script will take a list of file paths and either move and/or rename the files
as specified

=head1 OPTIONS

Database options

This is the database the objects will be written to

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the 
         standard port for mysql-g1kdcc.ebi.ac.uk

-move_list, tab file containing the to and from locations for the files you want 
            moved or renamed

-from, the file you wish moved or renamed
-to, the location/new name you wish to assign

-priority, the priority you wish to give to the LSF job running the actual move

-clobber, if you wish to clobber existing files

-run, must be switched on in order to store the actual objects

=head1 Examples

perl reseqtrack/scripts/file/move_archive_files.pl -dbhost myhost -dbuser myuser -dbpass ***** -dbport 4197 -dbname my_db -move move.list -run

=cut
