#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 delete_file);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use ReseqTrack::Tools::RunSamtools;
use File::Basename qw(fileparse);
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
my $samtools;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $delete_inputs;
my $reference;
my $reference_index;
my $directory_layout;
my $command;
my %options;
my $index_outputs;

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
  'samtools=s' => \$samtools,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'delete_inputs!' => \$delete_inputs,
  'reference=s' => \$reference,
  'reference_index=s' => \$reference_index,
  'directory_layout=s' => \$directory_layout,
  'command=s' => \$command,
  'options=s' => \%options,
  'index_outputs!' => \$index_outputs,
    );

my @allowed_cmds = qw(merge sort index fix_and_calmd calmd fixmate sam_to_bam);
throw("Don't recognise command $command. Acceptable commands are: @allowed_cmds")
  if (! grep {$command eq $_ } @allowed_cmds);

my @allowed_options = keys %{&ReseqTrack::Tools::RunSamtools::DEFAULT_OPTIONS};
foreach my $option (keys %options) {
  throw("Don't recognise option $option. Acceptable options are: @allowed_options")
    if (! grep {$option eq $_ } @allowed_options);
}

throw("Must specify an output directory") if (!$output_dir);
throw("Must specify an output type") if (!$type_output);
throw("Must specify an index type if index_outputs flag is used")
      if ($index_outputs && !$type_index);

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

if ($directory_layout) {
  my $rmia = $db->get_RunMetaInfoAdaptor;
  my $run_meta_info;
  if ($name =~ /[ESD]RR\d{6}/) {
    $run_meta_info = $rmia->fetch_by_run_id($&);
  }
  elsif ($name =~ /[ESD]RS\d{6}/) {
    my $rmi_list = $rmia->fetch_by_sample_id($&);
    $run_meta_info = $rmi_list->[0] if (@$rmi_list);
  }

  if ($run_meta_info) {
    $output_dir = create_directory_path($run_meta_info, $directory_layout, $output_dir);
  }
}

my $samtools_object = ReseqTrack::Tools::RunSamtools->new(
                    -input_files             => \@input_filepaths,
                    -program                 => $samtools,
                    -working_dir             => $output_dir,
                    -job_name                => $name,
                    -reference               => $reference,
                    -reference_index         => $reference_index,
                    -options                 => \%options,
                    );
$samtools_object->run($command);

my $output_bam_paths = ($command ne 'index')? $samtools_object->output_files : [];
my $index_file_paths = ($command eq 'index')? $samtools_object->output_files : [];
if ($index_outputs) {
    my $indexer = ReseqTrack::Tools::RunSamtools->new(
                    -input_files             => $output_bam_paths,
                    -program                 => $samtools,
                    );
    $indexer->run('index');
    $index_file_paths = $indexer->output_files;
}

if($store){
  my $host = get_host_object($host_name, $db);

  if (@$output_bam_paths) {
    foreach my $path (@$output_bam_paths) {
      throw("database already has file with name $path")
          if ($fa->fetch_by_name($path));
    }
    my $bams = create_objects_from_path_list($output_bam_paths, $type_output, $host);
    foreach my $bam (@$bams) {
      $bam->md5( run_md5($bam->name) );
    }
    my $collection = ReseqTrack::Collection->new(
        -name => $name, -type => $type_output,
        -others => $bams);
    $ca->store($collection);
  }

  if (@$index_file_paths) {
    my $fa = $db->get_FileAdaptor;
    my $bais = create_objects_from_pathlist($index_file_paths, $type_index, $host);
    foreach my $bai (@$bais) {
      $bai->md5( run_md5($bai->name) );
      $fa->store($bai);
    }
  }
}

if($delete_inputs){
  my @index_files;
  my @index_filepaths;
  INDEX_FILE:
  foreach my $filepath (@input_filepaths) {
    my $index_filepath = $filepath . '.bai';
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

reseqtrack/scripts/process/run_samtools.pl

=head1 SYNOPSIS

This script runs samtools to process sam / bam files.  It will do one of the following;
  
      merge: all sam / bam files in the collection are merged into a single sorted bam file
      sort: sorts each sam / bam file in a collection (does not merge)
      index: creates a bai file for each bam in a collection
      fix_and_calmd: Runs fixmate and calmd for each bam in a collection (does not merge)
      calmd: Runs calmd for each bam in a collection (does not merge)
      fixmate: Runs fixmate for each bam in a collection (does not merge)
      sam_to_bam: converts sam files to bam files (does not merge)

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

  -type_index, file type of the output bam index file
      (used if the -index_outputs flag is used or if command= 'index')

  -output_dir, base directory to hold files that do not need to be merged

  -samtools, the samtools executable.
  NB the Samtools class can guess the location of samtools if not specified.

  -host_name, default is '1000genomes.ebi.ac.uk', needed for storing output files

  -store, boolean flag, to store output files in the database.

  -delete_inputs, boolean flag, to delete the original input files (and index if it exists)
  and remove mark them as deleted in the database History table

  -reference, path to the reference fasta file.
  Used for calmd, and for converting sam to bam when the use_reference_index option is used

  -reference_index, path to the reference index file.
  Used for converting sam to bam when the use_reference_index option is used

  -directory_layout, specifies where the files will be located under output_dir.
      Tokens matching method names in RunMetaInfo will be substituted with that method's
      return value.

  -command, must be one of the following: 'merge', 'sort', 'index', 'fix_and_calmd', 'calmd', 'fixmate', 'sam_to_bam'
  Tells the RunSamtools object how to process the input files

  -options, for constructing the options hash passed to the RunSamtools object
  e.g. -options use_reference_index=1

  -index_outputs, flag to create an index file for any output bam file

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/process/run_samtools.pl $DB_OPTS -command calmd -index_outputs
    -type_input BAM -type_output MD_BAM -type_index MD_BAI
    -store -output_dir /path/to/base_dir/ -directory_layout population/sample_id/run_id
    -reference /path/to/reference -options input_sort_status=n

=cut

