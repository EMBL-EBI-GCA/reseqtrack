#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use File::Basename qw(dirname);
use ReseqTrack::Tools::RunValidateBam;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $name;
my $type_input;
my $type_output;
my $program;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $use_db_md5;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'name=s' => \$name,
  'type_input=s' => \$type_input,
  'type_output=s' => \$type_output,
  'program=s' => \$program,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'use_db_md5!' => \$use_db_md5,
    );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );
$db->dbc->disconnect_when_inactive(1);


my $ca = $db->get_CollectionAdaptor;
my $collection = $ca->fetch_by_name_and_type($name, $type_input);

throw("Failed to find a collection for ".$name." ".$type_input." from ".$dbname) 
    unless($collection);

my @input_filenames = map {$_->name} @{$collection->others};
my $working_dir = dirname($input_filenames[0]);

my $bam_validator = ReseqTrack::Tools::RunValidateBam->new(
                  -input_files             => \@input_filenames,
                  -working_dir             => $working_dir,
                  -program                 => $program,
                  -job_name                => $name,
                  );

if ($use_db_md5) {
  my $fa = $db->get_FileAdaptor;
  foreach my $file (@{$collection->others}) {
    $bam_validator->md5($file->name, $file->md5);
  }
}

$bam_validator->run;

$db->dbc->disconnect_when_inactive(0);
if($store){
  my $host = get_host_object($host_name, $db);
  my $fa = $db->get_FileAdaptor;

  my $bas_list = create_objects_from_path_list($bam_validator->output_files, $type_output, $host);

  foreach my $bas (@$bas_list) {
    $bas->md5( run_md5($bas->name) );
    $fa->store($bas);
  }

}


=pod

=head1 NAME

reseqtrack/scripts/process/run_squeeze_bam.pl

=head1 SYNOPSIS

This script is used reduce the file size of a bam using the squeeze function of bamUtil

File size is reduced by optionally:
  dropping OQ fields, using the -rm_OQ_tags command line option
  dropping any other field, using the rm_tag_type command line option
  removing duplicates, using the rm_dups command line option

The input bam file is taken from a collection in the database. The output bam will be written to the database.
The input bam file can be deleted, along with its index file, and this will be recorded in the History table of the database.

=head1 OPTIONS

  database options:

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on

  other options:

  -rm_OQ_fields, boolean flag, remove OQ_fields from the bam

  -rm_dups, boolean flag, remove reads marked as library duplicates

  -rm_tag_type, string, e.g. 'XM:i', can be specified multiple times, remove these fields from the bam

  -name, name of a collection containing the bam file
  If name is a run_id (or contains a run_id), the run_meta_info table will be used to get some info

  -type_input, type of the collection containing the input bam file

  -type_output, collection type and file type when storing output bam file in the database

  -output_dir, base directory to hold files that do not need to be merged

  -program, path to the bamUtil 'bam' executable.

  -store, boolean flag, to store output bam file in the database.

  -disable_md5, boolean flag, files written to the database will not have an md5

  -host_name, default is '1000genomes.ebi.ac.uk', needed for storing output bam file

  -delete_input, boolean flag, to delete the original bam file (and index if it exists)
  and remove mark them as deleted in the database History table

  -directory_layout, specifies where the files will be located under output_dir.
      Tokens matching method names in RunMetaInfo will be substituted with that method's
      return value.

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/process/run_squeeze_bam.pl $DB_OPTS -name SRR000001
    -rm_OQ_fields -rm_tag_type XM:i -rm_tag_type XG:i
    -type_input BAM -type_output SQUEEZED_BAM -output_dir /path/to/base_dir/
    -store -delete_input -directory_layout population/sample_id/run_id
    -program /path/to/executable

=cut

