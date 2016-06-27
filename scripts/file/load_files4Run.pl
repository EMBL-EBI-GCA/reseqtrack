#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use ReseqTrack::Tools::Loader::File;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::DBSQL::DBAdaptor;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $dir;
my $list_file;
my $type;
my $host_name = "1000genomes.ebi.ac.uk";
my $run;
my $die_for_problems = 1;
my $update_existing;
my $store_new;
my $do_md5=1;
my $md5_program = 'md5sum';
my $help;
my $assign_types = 1;

&GetOptions(
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'list_file=s' => \$list_file,
  'type|file_type=s' => \$type,
  'host=s' => \$host_name,
  'run!' => \$run,
  'die_for_problems!' => \$die_for_problems,
  'update_existing!' => \$update_existing,
  'store_new!' => \$store_new,
  'md5_program=s' => \$md5_program,
  'help!' => \$help,
    );

if($help){
    usage();
}

usage() if !$list_file || !$type || !$dbhost || !$dbname || !$dbuser || !$dbpass || !$host_name;

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
                        -host   => $dbhost,
                        -user   => $dbuser,
                        -port   => $dbport,
                        -dbname => $dbname,
                        -pass   => $dbpass,
 );

throw( "No DB connection established" ) if ( !$db );

my $file_paths=get_filepaths_per_run($list_file);

foreach my $run_id (keys %{$file_paths}) {
    my @files=@{$file_paths->{$run_id}};

    my $loader = ReseqTrack::Tools::Loader::File->new(
	-dir       => $dir,
	-file      => \@files,
	-type      => $type,
	-hostname      => $host_name,
	-dbhost => $dbhost,
	-dbname => $dbname,
	-dbuser  => $dbuser,
	-dbpass  => $dbpass,
	-dbport  => $dbport,
	-do_md5  => $do_md5,
	-update_existing=>$update_existing,
	-store_new=>$store_new,
	);

    $loader->process_input();
    $loader->create_objects();
    $loader->sanity_check_objects();
    $loader->load_objects() if($run);

    print STDOUT "[INFO] Creating ReseqTrack::Collection with name $run_id and type $type\n";

    my $collection = ReseqTrack::Collection->new(
        -name => $run_id,
        -type => $type,
        -others => $loader->objects,
        -table_name =>'file',
	);
    $db->get_CollectionAdaptor->store($collection) if $run;
}

sub usage{
    exec('perldoc', $0);
    exit(0);
}

sub get_filepaths_per_run {
    my $file_list=shift;

    my %files;

    open FH,"<$file_list" or throw("Can't open $file_list: $!");

    while(<FH>) {
        chomp;
        my ($run_source_id,$file)=split/\t/;
        throw("Error reading $file_list. I need two tab-separated columns: run_source_id\\tfilepath") unless ($run_source_id && $file);
        push @{$files{$run_source_id}},$file;
    }

    close FH;

    return \%files;
}

=pod

=head1 NAME

ReseqTrack/scripts/files/load_files4Run.pl

=head1 SYNOPSIS

This script will take a file with a list of file paths for each of the run ids and will load
the file table of a ReseqTrack database with the details and will create a new entry in the Collection for each run.

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

-list_file, a file containing a list of run ids and files paths in the format: run_id\tfile_path

-md5_file, a list of md5sums for each file in the format: md5checkum filepath.

Other options

-type, This should be a string which will be associated with all files and with the new collections. For example FASTQ if you are loading fastq files. There are no restrictions on what this is other than it should be shorter than 50characters, convention normally has each type in 
upper case and it is good if it is in someway informative about the files loaded. If you specify the -assign_types option you don't have to use this

Note: md5 checksums will be calculated for each file to be loaded.

-host, This is the name of the host which the filesystem is visible to so 1000genomes.ebi.ac.uk for ebi files.

-run, This is a binary option which needs no additional argument when specified the 
adaptor store method is run. By default this is off.

-update_existing, This means if the store method finds a file with the same name it
will use the update method rather than the store method

-store_new, This means if the store method finds a file with the same name it will
still store a brand new line unless the md5s are the same

-md5_program, This is a string pointing to the md5sum program

-help, This makes the script print out its options

=head1 Examples

Run it by passing a md5.txt file containing the md5sum per file:

perl load_files4Run.pl -list_file filesPerRun.txt -dbhost mysql-blueprint -dbuser g1krw -dbpass XXXXX -dbport 4360 -dbname test_db -type FASTQ -host 1000genomes.ebi.ac.uk -run

Where filesPerRun.txt could be (in this case 2 paired-end FASTQ files for the same run id):

    SRXXX    /path/to/file/SRRXXX_1.fastq
    SRXXX    /path/to/file/SRRXXX_2.fastq

=cut

