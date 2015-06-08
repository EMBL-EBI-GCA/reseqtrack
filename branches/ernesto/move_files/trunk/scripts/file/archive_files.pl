#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;
use ReseqTrack::Tools::Loader;
use ReseqTrack::Tools::Loader::Archive;
use ReseqTrack::Tools::GeneralUtils;

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
my $no_lock;
my $help;


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
	    'from_db!'        => \$from_db,
	    'path_like:s'     => \$path_like,
	    'archive_sleep=s' => \$sleep, 
	    'no_lock!'        => \$no_lock,
	    'help!' => \$help,
	   );

if($help){
  useage();
}

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




=pod

=head1 NAME

ReseqTrack/scripts/files/archive_files.pl

=head1 SYNOPSIS

This script will take a list of file paths from a variety of options and load
the archive table with them and the given action so they can be moved onto or off
the fire archive disk.

=head1 OPTIONS

Database options

This is the database the objects will be written to

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the standard
port for mysql-g1kdcc.ebi.ac.uk

File list and other options

-file, the path to a single file, this can appear on the commandline multiple times.

-list_file, a file containing a list of files paths, one to each line,

-dir, If no other options are specified all the files are read from this directory
and the directory tree below it

-path_like, use the the given string to fetch_all_by_path_like from the file table

-type, The type of file to retrieve from the database

-from_db, This indicates that the paths must be retrieved from the database either using the specifed path, or the type or both

-run, This is a binary option which needs no additional argument when specified the 
adaptor store method is run. By default this is off.

-descend, If a directory is passed in this means the directory tree is descended and
all files below that point are stored. By default this option in on but you can use
the -nodescend option to prevent it happening

-action, which archive action to do, this can be worked out automatically from the 
path 

-priority, the priority to give the LSF job which actually runs the archiving process

-skip, this stops the script for continually looping over the archive table to clean it up

-max_number_of_archives, this specifies the maximum number of archive objects that can be in the system at once, by default this is 1000, this is to stop the system from getting to bogged down

-archive_sleep, this is the amount of time to sleep for between cleanup cycles, by default it is 240 seconds

-no_lock, this switches off the checking of the database meta table for a locking string

-help, This makes the script print out its options

=head1 Examples

perl reseqtrack/scripts/file/archive_files.pl -dbhost my_host -dbuser username -dbpass ****** -dbport 4197 -dbname my_database -action archive -skip -run -priority 90 -from_db -type FILTERED_FASTQ

This would mark for archive all FILTERED_FASTQ files in the staging area

=cut
