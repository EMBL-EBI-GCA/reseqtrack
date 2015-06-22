#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_object_from_path);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 );
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunTransposeBam;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $name;
my $type_input_string;
my $type_input_collection;
my $type_output;
my $output_dir;
my $program;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $disable_md5;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'name=s' => \$name,
  'type_input_string=s' => \$type_input_string,
  'type_input_collection=s' => \$type_input_collection,
  'type_output=s' => \$type_output,
  'output_dir=s' => \$output_dir,
  'program=s' => \$program,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'disable_md5!' => \$disable_md5,
    );

throw("Must specify an output directory") if (!$output_dir);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );
$db->dbc->disconnect_when_inactive(1);


my $is_a = $db->get_InputStringAdaptor;
my $input_string = $is_a->fetch_by_name_and_type($name, $type_input_string);
throw("Failed to find an input strint for $name and $type_input_string") if (!$input_string);

my $ca = $db->get_CollectionAdaptor;
my ($col_name, $region) = split(/\./, $name);
throw("did not find a collection name from input string $name") if !$col_name;
my $collection = $ca->fetch_by_name_and_type($col_name, $type_input_collection);
throw("did not find a collection for name $col_name and type $type_input_collection") if !$collection;

my $files = $collection->others;
throw("no files for collection") if !@$files;
my @filenames = map{$_->name} @$files;

$output_dir .= "/$col_name";
$output_dir =~ s{//+}{/}g;

my $job_name = $col_name . '.' . $region;
$job_name =~ s/:/./;


my $bam_transposer = ReseqTrack::Tools::RunTransposeBam->new(
                  -input_files             => \@filenames,
                  -working_dir             => $output_dir,
                  -program                 => $program,
                  -job_name                => $job_name,
                  -region                  => $region,
                  );

$bam_transposer->run;

if($store){
  my $host = get_host_object($host_name, $db);
  my $fa = $db->get_FileAdaptor;

  my $bam = create_object_from_path($bam_transposer->output_files->[0], $type_output, $host);
  if (! $disable_md5) {
    $bam->md5( run_md5($bam->name) );
  }
  $fa->store($bam);
}


=pod

=head1 NAME

reseqtrack/scripts/process/run_transpose_bam.pl

=head1 SYNOPSIS

This script is used to merge sample bams and simultaneously split by region

The input string is name after the collection of bams and the region, separated by a '.'
  e.g. 'GBR_bwa.chr5:1000-1000' or 'GBR_bwa.chr5'
  where GBR_bwa is the name of a collection of bam files

=head1 OPTIONS

  database options:

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on

  other options:

  -name, an input string loaded into the input string table, as described above

  -type_input_string, type of the input_string in the database

  -type_input_collection, type of the bam file collection in the database

  -type_output, type used for storing output files

  -output_dir, base directory for output files

  -program, path to the transpose_bam executable in reseqtrack/c_code

  -store, boolean flag, to store output bam file in the database.

  -disable_md5, boolean flag, files written to the database will not have an md5

  -host_name, default is '1000genomes.ebi.ac.uk', needed for storing output bam file

=cut

