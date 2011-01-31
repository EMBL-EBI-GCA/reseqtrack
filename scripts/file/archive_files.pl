#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use ReseqTrack::Tools::Loader;
use ReseqTrack::Tools::Loader::Archive;
use Data::Dumper;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $dir;
my @files;
my $list_file;
my $remote;
my $action;
my $run;
my $descend= 1;
my $debug;
my $from_db;
my $type;
my $path_like;
my $action_string = 'archive';
my $action_location_name;
my $sleep = 240;
my $skip_cleanup = 0;
my $verbose = 0;
my $priority = 50;
my $max_number = 1000;
my $lines_check = 1;
my $no_lock;



&GetOptions(
	    'dbhost=s'               => \$dbhost,
	    'dbname=s'               => \$dbname,
	    'dbuser=s'               => \$dbuser,
	    'dbpass=s'               => \$dbpass,
	    'dbport=s'               => \$dbport,
	    'dir=s'                  => \$dir,
	    'file=s@'                => \@files,
	    'list_file|file_list=s'  => \$list_file,
            'type|file_type=s'       => \$type,
	    'run!'                   => \$run,
	    'verbose!'               => \$verbose,
	    'descend!'               => \$descend,
	    'debug!'                 => \$debug,
	    'action=s'               => \$action_string,
	    'priority=s'             =>\$priority,
	    'skip!'                  =>\$skip_cleanup,
	    'max_number_of_archives=s' => \$max_number,
	    'lines_check!'    => \$lines_check,	    
	    'from_db!'        => \$from_db,
	    'type:s'          => \$type,
	    'path_like:s'     => \$path_like,
	    'archive_sleep=s' => \$sleep, 
	    'no_lock!'        => \$no_lock,
	   );

my $archiver = ReseqTrack::Tools::Loader::Archive->new(
						       -dir       => $dir,
						       -file      => \@files,
						       -list_file => $list_file,
						       -type      => $type,
						       -descend   => $descend,
						       -dbhost    => $dbhost,
						       -dbname    => $dbname,
						       -dbuser    => $dbuser,
						       -dbpass    => $dbpass,
						       -dbport    => $dbport,
						       -debug     => $debug,
						       -action    => $action_string,
						       -verbose   => $verbose,
						       -priority  => $priority,
						       -max_number=> $max_number,
						       -from_db   => $from_db,
						       -path_like => $path_like,
						       -archive_sleep => $sleep,
						       -no_lock   =>$no_lock,
						      );
$archiver->process_input();
#$archiver->cleanup_archive_table($verbose);
$archiver->sanity_check_objects();
$archiver->archive_objects() if $run;


my $max_tries = 10;
my $tries     =  0;
if (!$skip_cleanup && $run) {
  print "Attempting cleanup of archive table in 60 seconds\n";
  sleep (60);
  my $clean_archive_table = 1;
  while ($clean_archive_table) {
    $tries++;
    my $obs_remaining = $archiver->cleanup_archive_table($verbose);
    if ($obs_remaining) {
     
      print "Found $obs_remaining. Waiting for $sleep\n";
      sleep($sleep);
    } else {
      $clean_archive_table = 0;
    }
  }
}




