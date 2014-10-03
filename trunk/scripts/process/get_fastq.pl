#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::FileUtils qw(create_object_from_path create_history);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::ERAUtils qw(get_erapro_conn);
use ReseqTrack::Tools::MetaDataUtils qw( create_directory_path );

$| = 1;

my $run_id;
my $source_root_location;
my $host_name = '1000genomes.ebi.ac.uk';
my $file_type;
my $clobber;
my $dbhost;
my $dbname;
my $dbuser;
my $dbport;
my $dbpass;
my $output_dir;
my $load = 0;
my $era_dbuser;
my $era_dbpass;
my $era_dbname;
my $help = 0;
my $directory_layout = 'sample_alias/archive_sequence';
my $module = 'ReseqTrack::Tools::GetFastq';
my %module_options;

&GetOptions( 
  'run_id=s' => \$run_id,
  'output_dir=s' => \$output_dir,
  'host_name=s' => \$host_name,
  'type|file_type=s' => \$file_type,
  'clobber!' => \$clobber,
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'load!' => \$load,
  'era_dbuser=s' =>\$era_dbuser,
  'era_dbpass=s' => \$era_dbpass,
  'era_dbname=s' => \$era_dbname,
  'help!' => \$help,
  'directory_layout=s' => \$directory_layout,
  'module=s' => \$module,
  'root_location=s' => \$source_root_location,
  'module_option=s' => \%module_options,
    );

if($help){
  usage();
}

eval "require $module" or throw "cannot load module $module $@";

throw("need a run_id") if !$run_id;
throw("need a file type to load files") if ($load && !$file_type);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $era_db = get_erapro_conn($era_dbuser, $era_dbpass, $era_dbname);

my $run = $db->get_RunAdaptor->fetch_by_source_id($run_id);

throw("Failed to find a run object for ".$run_id." from ".$dbname)
    unless($run);

if($run->instrument_platform && $run->instrument_platform eq 'COMPLETE_GENOMICS'){
  exit(0);
}

if ($directory_layout) {
  $output_dir = create_directory_path($run, $directory_layout, $output_dir);
}

$db->dbc->disconnect_when_inactive(1);

my %constructor_hash;
while (my ($key, $value) = each %module_options) {
  $constructor_hash{'-'.$key} = $value;
}
$constructor_hash{-output_dir} = $output_dir;
$constructor_hash{-run_info} = $run;
$constructor_hash{-source_root_dir} = $source_root_location;
$constructor_hash{-clobber} = $clobber;
$constructor_hash{-db} = $era_db;

my $fastq_getter = $module->new (%constructor_hash);
my $num_files = $fastq_getter->run;

my $host = get_host_object($host_name, $db);
my $fastq_paths = $fastq_getter->output_files;
my @file_objects;
if($load && $num_files){
  $db->dbc->disconnect_when_inactive(0);
  my $fa = $db->get_FileAdaptor;
  foreach my $path (@$fastq_paths) {
    my $file = create_object_from_path($path, $file_type, $host);
    my $md5 = $fastq_getter->md5_hash->{$path};
    warning("Don't have md5 for $path") if !$md5;
    $file->md5($md5);

    my $stored_files = $fa->fetch_by_filename($file->filename);
    if (@$stored_files) {
      if (@$stored_files > 1) {
        throw("Not sure what to do, seem to have ".@$stored_files." associated with ".  $file->filename)
      }
      my $history = create_history($file, $stored_files->[0]);
      $file->dbID($stored_files->[0]->dbID);
      $file->history($history);
    }
    push(@file_objects, $file);
  }

  my $collection =  ReseqTrack::Collection->new(
    -name => $run_id,
    -others => \@file_objects,
    -type => $file_type,
      );

  my $ca = $db->get_CollectionAdaptor;
  $ca->store($collection);

  my $stored_collection = $ca->fetch_by_name_and_type($run_id, $file_type);
  foreach my $stored_file(@{$stored_collection->others}){
    throw($stored_file->name." doesn't exist") unless(-e $stored_file->name);
  }
}




sub usage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/process/get_fastq.pl

=head1 SYNOPSIS

This script fetches fastq files from the ERA ftp servers.

=head1 OPTIONS

Database options

This is the database the objects will be written to

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the standard
port for mysql-g1kdcc.ebi.ac.uk

These are for the ERA database

-era_dbuser, the name of the era database user
-era_dbpass, the password for the era database

Controls method of getting fastq:

-module, default is ReseqTrack::Tools::GetFastq.  Can use any module that inherits from this class.
-module_option, used only for child classes of GetFastq, e.g. -option password=mypassword

Input options

-run_id, this is the run id for the run you wish to fetch the fastq files for
-output_dir, this is the root path where the files should be written. Meta information associated with the run id will be used to complete the path in the style root_path/sample_id/archive_fastq
-type/file_type, this is the file/collection type to use when loading the data
-host_name, this is the name of the host object to associated with the files (default 1000genomes.ebi.ac.uk)
-clobber, this specifies whether existing files should be copied over or not
-load, this species if the created file/collection objects should be loaded into
the database
-directory_layout, this specifies where the files will be located under output_dir. Tokens matching method names in RunMetaInfo will be substituted with that methods return value. Default value is sample_name/archive_sequence
-root_location, the root directory for era files, default is /nfs/era-pub

=head1 Examples

perl ReseqTrack/scripts/process/get_fastq.pl -output scratch/staging_area/sequence_staging/ -clobber -dbhost host -dbuser user -dbpass **** -dbport 4197 -dbname staging_db -host 1000genomes.ebi.ac.uk -file_type ARCHIVE_FASTQ -load -run_id SRR031624 -era_dbuser username -era_dbpass *****

=cut

