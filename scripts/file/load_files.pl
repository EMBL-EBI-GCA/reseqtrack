#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use ReseqTrack::Tools::Loader;
use ReseqTrack::Tools::Loader::File;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list create_history assign_type check_type);
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::IndexUtils qw(get_md5hash_from_sequence_index get_md5hash_from_alignment_index);
use ReseqTrack::DBSQL::DBAdaptor;

use ReseqTrack::File;
use ReseqTrack::Host;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $dir;
my @files;
my $list_file;
my $md5_file;
my $type;
my $host_name;
my $run;
my $verbose;
my $descend = 1;
my $die_for_problems = 1;
my $update_existing;
my $store_new;
my $check_md5;
my $md5_program = 'md5sum';
my $help;
my $assign_types = 1;
my $check_types = 1;
my $debug = 0;
&GetOptions(
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'dir=s' => \$dir,
  'file=s@' => \@files,
  'list_file=s' => \$list_file,
  'md5_file=s' => \$md5_file,
  'type|file_type=s' => \$type,
  'host=s' => \$host_name,
  'run!' => \$run,
  'verbose!' => \$verbose,
  'descend!' => \$descend,
  'die_for_problems!' => \$die_for_problems,
  'update_existing!' => \$update_existing,
  'store_new!' => \$store_new,
  'check_md5!' => \$check_md5,
  'md5_program=s' => \$md5_program,
  'help!' => \$help,
  'assign_types!' => \$assign_types,
  'check_types!' => \$check_types,
  'debug' =>\$debug, 
   );

if($help){
  useage();
}

if(!$type && !$assign_types){
  throw("Must give load_files.pl a file type with -type");
}

if(!$host_name){
  throw("Must give load_files a host name with -host");
}

if(!$dbhost || !$dbname || !$dbuser){
  throw("Must provide database connection details with -dbhost -dbuser -dbpass ".
        "-dbport -dbname");
}


my $loader = ReseqTrack::Tools::Loader::FileInherit->new(
 -dir       => $dir,
 -file      => \@files,
 -list_file => $list_file,
 -type      => $type,
 -host      => $host_name,
  -md5_file  => $md5_file,
 -hostname      => $host_name,
 -dbhost => $dbhost,
 -dbname => $dbname,
 -dbuser  => $dbuser,
 -dbpass  => $dbpass,
 -dbport  => $dbport,
 -debug => $debug,
  -assign_types => $assign_types,
  -check_types=>$check_types,
  -check_md5  => $check_md5,
  -update_existing=>$update_existing,
  -store_new=>$store_new,
  -debug=>$debug,
);

$loader->process_input();
$loader->create_files();
$loader->load_files($run);





sub useage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/files/load_files.pl

=head1 SYNOPSIS

This script will take a list of file paths from a variety of options and load
the file table of a ReseqTrack database with the details

=head1 OPTIONS

Database options

This is the database the objects will be written to

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the standard
port for mysql-g1kdcc.ebi.ac.uk

File list options

-file, the path to a single file, this can appear on the commandline multiple times.

-list_file, a file containing a list of files paths, one to each line,

-md5_file, a file of output from md5sum in the format md5checkum filepath.

-alignment_index_file, a file in the format of the ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README.alignment.index, both the bam and the bai files will be stored

-sequence_index_file, a file in the format ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README.sequence.index.

-dir, If no other options are specified all the files are read from this directory
and the directory tree below it and stores. If you only want to store files in
that directory use -nodescend option.

-use_dir, If this is specified along with -dir string when any of the other file lists
are passed in the -dir string is associated with each path string given.

Other options

-type, This should be a string which will be associated with all files for example FASTQ if you are loading fastq files. There are no restrictions on what this is other 
than it should be shorter than 50characters, convention normally has each type in 
upper case and it is good if it is in someway informative about the files loaded. If you specify the -assign_types option you don't have to use this

-host, This is the name of the host which the filesystem is visible to so 1000genomes.ebi.ac.uk for ebi files.

-run, This is a binary option which needs no additional argument when specified the 
adaptor store method is run. By default this is off.

-verbose, This switchs on some print statements about the process. By default this is
off

-descend, If a directory is passed in this means the directory tree is descended and
all files below that point are stored. By default this option in on but you can use
the -nodescend option to prevent it happening

-die_for_problems, This is a binary option which alters how the script behaves when
it seems problems in the given file list. If this is switched on it prints out
all problems then throws and exception if not it just prints out the problems. By default this is off

-update_existing, This means if the store method finds a file with the same name it
will use the update method rather than the store method

-store_new, This means if the store method finds a file with the same name it will
still store a brand new line unless the md5s are the same

-check_md5, This makes the script run md5 checksums on each file

-md5_program, This is a string pointing to the md5sum program

-assign_types, This uses the ReseqTrack::Tools::FileUtils assign_type method to 
specify file types on the basis of the file path, this is on by default but it can be switched off by using -noassign_types

-check_types, This checks the type of each file object against an acceptable list, this is on by default but it can be switched off with -nocheck_types

-help, This makes the script print out its options

=head1 Examples

To load a single file

perl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname lec_test_track -file_type MISC -noassign_type -host 1000genomes.ebi.ac.uk -file /path/to/file -run

To load from an index file

perl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname lec_test_track  -host 1000genomes.ebi.ac.uk -sequence_index_file /path/to/sequence.index -run -use_dir -dir /root/to/sequence/files

To load from a directory

perl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname lec_test_track -host 1000genomes.ebi.ac.uk -dir /path/to/dir -run


=cut
