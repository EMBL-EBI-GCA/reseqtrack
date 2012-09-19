#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_object_from_path);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 delete_file);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use File::Basename qw(fileparse);
use ReseqTrack::Tools::RunBamSqueeze;
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
my $output_dir;
my $program;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $disable_md5;
my $delete_input;
my $directory_layout;
my $rm_OQ_fields;
my $rm_dups;
my @rm_tag_types;
my $run_id_regex = '[ESD]RR\d{6}';
my $sample_id_regex = '[ESD]RS\d{6}';

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'name=s' => \$name,
  'type_input=s' => \$type_input,
  'type_output=s' => \$type_output,
  'output_dir=s' => \$output_dir,
  'program=s' => \$program,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'disable_md5!' => \$disable_md5,
  'delete_input!' => \$delete_input,
  'directory_layout=s' => \$directory_layout,
  'rm_OQ_fields!' => \$rm_OQ_fields,
  'rm_dups!' => \$rm_dups,
  'rm_tag_type=s' => \@rm_tag_types,
  'run_id_regex=s' => \$run_id_regex,
    );

throw("Must specify an output directory") if (!$output_dir);
throw("Must specify one or more of -rm_OQ_fields, -rm_dups, -rm_tag_type")
    if (! grep {$_} ($rm_OQ_fields, $rm_dups, scalar @rm_tag_types));

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

my $others = $collection->others;
throw("Expecting one file; have ".@$others."\n")
    if (@$others != 1);
my $input_file = $others->[0];


if ($directory_layout) {
  my $rmia = $db->get_RunMetaInfoAdaptor;
  my $run_meta_info;
  if ($name =~ /$run_id_regex/) {
    $run_meta_info = $rmia->fetch_by_run_id($&);
  }
  elsif ($name =~ /$sample_id_regex/) {
    my $rmi_list = $rmia->fetch_by_sample_id($&);
    $run_meta_info = $rmi_list->[0] if (@$rmi_list);
  }

  if ($run_meta_info) {
    $output_dir = create_directory_path($run_meta_info, $directory_layout, $output_dir);
  }
}


my $bam_squeezer = ReseqTrack::Tools::RunBamSqueeze->new(
                  -input_files             => $input_file->name,
                  -working_dir             => $output_dir,
                  -program                 => $program,
                  -job_name                => $name,
                  -rm_tag_types            => \@rm_tag_types,
                  );

$bam_squeezer->options('keepOQ', $rm_OQ_fields ? 0 : 1);
$bam_squeezer->options('keepDups', $rm_dups ? 0 : 1);

$bam_squeezer->run;

if($store){
  my $host = get_host_object($host_name, $db);

  my $bam = create_object_from_path($bam_squeezer->output_files->[0], $type_output, $host);
  if (! $disable_md5) {
    $bam->md5( run_md5($bam->name) );
  }
  my $collection = ReseqTrack::Collection->new(
      -name => $name, -type => $type_output,
      -others => $bam);
  $ca->store($collection);

}

if($delete_input){
  my @db_files = ($input_file);
  my @delete_filepaths = ($input_file->name);
  my $ha = $db->get_HistoryAdaptor;
  my $fa = $db->get_FileAdaptor;

  my $index_filepath = $input_file->name . '.bai';
  if (-f $index_filepath) {
    push(@delete_filepaths, $index_filepath);
    my $index_file = $fa->fetch_by_name($index_filepath);
    if ($index_file) {
      push(@db_files, $index_file);
    }
  }

  foreach my $file (@db_files) {
    my $history = ReseqTrack::History->new(
        -other_id => $file->dbID, -table_name => 'file',
        -comment => "deleted by $0");
    $ha->store($history);
  }
  foreach my $filepath (@delete_filepaths) {
    delete_file($filepath);
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

  -run_id_regex, used to get run meta info.  Default is '[ESD]RR\d{6}'
  -study_id_regex, used to get run meta info.  Default is '[ESD]RS\d{6}'

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/process/run_squeeze_bam.pl $DB_OPTS -name SRR000001
    -rm_OQ_fields -rm_tag_type XM:i -rm_tag_type XG:i
    -type_input BAM -type_output SQUEEZED_BAM -output_dir /path/to/base_dir/
    -store -delete_input -directory_layout population/sample_id/run_id
    -program /path/to/executable

=cut

