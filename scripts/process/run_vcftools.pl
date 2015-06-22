#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 delete_file);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunVcfTools;
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
my $type_index;
my $output_dir;
my $vcftools_dir;
my $reference_index;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $disable_md5;
my $delete_inputs;
my $command;
my %options;
my $create_index;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'name=s' => \$name,
  'type_input=s' => \$type_input,
  'type_output=s' => \$type_output,
  'type_index=s' => \$type_index,
  'output_dir=s' => \$output_dir,
  'vcftools_dir=s' => \$vcftools_dir,
  'reference_index=s' => \$reference_index,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'disable_md5!' => \$disable_md5,
  'delete_inputs!' => \$delete_inputs,
  'command=s' => \$command,
  'options=s' => \%options,
  'create_index!' => \$create_index,
);

my @allowed_cmds = ReseqTrack::Tools::RunVcfTools->get_valid_commands;
throw("Don't recognise command $command. Acceptable commands are: @allowed_cmds")
  if (! grep {$command eq $_ } @allowed_cmds);

my @allowed_options = keys %{&ReseqTrack::Tools::RunVcfTools::DEFAULT_OPTIONS};
foreach my $option (keys %options) {
  throw("Don't recognise option $option. Acceptable options are: @allowed_options")
    if (! grep {$option eq $_ } @allowed_options);
}

throw("Must specify an output directory") if (!$output_dir);
throw("Must specify an output type") if (!$type_output);
throw("Must specify an index type if create_index flag is used")
      if ($create_index && !$type_index);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );
$db->dbc->disconnect_when_inactive(1);
my $ca = $db->get_CollectionAdaptor;
my $fa = $db->get_FileAdaptor;

my $collection = $ca->fetch_by_name_and_type($name, $type_input);
throw("Failed to find a collection for ".$name." ".$type_input." from ".$dbname) 
    unless($collection);

my $input_files = $collection->others;
my @input_filepaths = map {$_->{'name'}} @$input_files;

my $vcftools_object = ReseqTrack::Tools::RunVcfTools->new(
                    -input_files             => \@input_filepaths,
                    -working_dir             => $output_dir,
                    -job_name                => $name,
                    -options                 => \%options,
                    -vcftools_dir            => $vcftools_dir,
                    -create_index            => $create_index,
                    -reference_index         => $reference_index,
                    );
$vcftools_object->run($command);

if($store){
  my $host = get_host_object($host_name, $db);

  my @file_paths = @{$vcftools_object->output_vcf_files};
    
  if (@file_paths) {
    foreach my $path (@file_paths) {
      throw("database already has file with name $path")
          if ($fa->fetch_by_name($path));
    }
    my $files = create_objects_from_path_list(\@file_paths, $type_output, $host);
    foreach my $file (@$files) {
      if (! $disable_md5) {
        $file->md5( run_md5($file->name) );
      }
      $fa->store($file);
    }
  }

  my $fa = $db->get_FileAdaptor;
  my $tbi_paths = $vcftools_object->output_tbi_files;
  if (@$tbi_paths) {
    my $tbis = create_objects_from_path_list($tbi_paths, $type_index, $host);
    foreach my $tbi (@$tbis) {
      if (! $disable_md5) {
        $tbi->md5( run_md5($tbi->name) );
      }
      $fa->store($tbi);
    }
  }
}

if($delete_inputs){
  my @index_files;
  my @index_filepaths;
  INDEX_FILE:
  foreach my $filepath (@input_filepaths) {
    my $index_filepath = $filepath . '.tbi';
    next INDEX_FILE if (! -f $index_filepath);
    push(@index_filepaths, $index_filepath);
    my $index_file = $fa->fetch_by_name($index_filepath);
    if ($index_file) {
      push(@index_files, $index_file);
    }
  }

  my $ha = $db->get_HistoryAdaptor;
  foreach my $file (@$input_files, @index_files) {
    my $history = ReseqTrack::History->new(
        -other_id => $file->dbID, -table_name => 'file',
        -comment => "deleted by $0");
    $ha->store($history);
  }
  foreach my $file (@input_filepaths, @index_filepaths) {
    delete_file($file);
  }
}

=pod

=head1 NAME

reseqtrack/scripts/process/run_vcftools.pl

=head1 SYNOPSIS

This script runs vcftools to process vcf files.  It will do one of the following;
  
      merge: all vcf files in the collection are merged into a single vcf file
      concat: all vcf files in the collection are concatenated together

The input files are taken from a collection in the database. The output files will be written to the database.
The input files can be deleted, along with any index files, and this will be recorded in the History table of the database.


=head1 OPTIONS

  database options:

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on

  other options:

  -name, name of the collection of input files
  If name is a run_id / sample_id (or contains a run_id / sample_id), the run_meta_info table will be used to get some info

  -type_input, type of the collection of input files

  -type_output, collection type and file type when storing output files in the database

  -type_index, file type of the output bam index file (used if the -create_index flag is used)

  -output_dir, base directory to hold files that do not need to be merged

  -vcftools_dir, path to the directory containing the vcftools executables

  -host_name, default is '1000genomes.ebi.ac.uk', needed for storing output files

  -reference_index, path to a .fai file. Only used if concatenating / sorting over multiple chromosomes

  -store, boolean flag, to store output files in the database.

  -disable_md5, boolean flag, files written to the database will not have an md5

  -delete_inputs, boolean flag, to delete the original input files (and index if it exists)
  and mark them as deleted in the database History table

  -command, must be one of the following: 'merge', 'concat'
  Tells the RunVcfTools object how to process the input files

  -options, for constructing the options hash passed to the VcfTools object
  e.g. -options remove_duplicates=1

  -create_index, flag to create a tabix index file for any output vcf file

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/process/run_vcftools.pl $DB_OPTS -command merge -create_index -name mycalls
    -type_input VCF -type_output MERGED_VCF -type_index MERGED_TBI
    -store -output_dir /path/to/output_dir/
    -options remove_duplicates => 1

=cut

