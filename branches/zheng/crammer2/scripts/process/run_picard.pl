#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 delete_file);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use ReseqTrack::Tools::RunPicard;
use ReseqTrack::Tools::StatisticsUtils;
use ReseqTrack::Tools::Loader::File;
use File::Basename qw(fileparse);
use Getopt::Long;
use Data::Dumper;
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
my $picard_dir;
my $java_exe;
my $jvm_options;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $disable_md5;
my $delete_inputs;
my $directory_layout;
my $run_id_regex = '[ESD]RR\d{6}';
my $sample_id_regex = '[ESD]RS\d{6}';
my $command;
my %options;
my $create_index;
my $store_stats;
my $metrics_file_type;

&GetOptions(
  'dbhost=s'           => \$dbhost,
  'dbname=s'           => \$dbname,
  'dbuser=s'           => \$dbuser,
  'dbpass=s'           => \$dbpass,
  'dbport=s'           => \$dbport,
  'name=s'             => \$name,
  'type_input=s'       => \$type_input,
  'type_output=s'      => \$type_output,
  'type_index=s'       => \$type_index,
  'output_dir=s'       => \$output_dir,
  'picard_dir=s'       => \$picard_dir,
  'java_exe=s'         => \$java_exe,
  'jvm_options=s'      => \$jvm_options,
  'host_name=s'        => \$host_name,
  'store!'             => \$store,
  'disable_md5!'       => \$disable_md5,
  'delete_inputs!'     => \$delete_inputs,
  'directory_layout=s' => \$directory_layout,
  'command=s' => \$command,
  'options=s' => \%options,
  'run_id_regex=s' => \$run_id_regex,
  'sample_id_regex=s' => \$sample_id_regex,
  'create_index!' => \$create_index,
  'store_stats!' => \$store_stats,
  'metrics_type=s' => \$metrics_file_type,
);

my @allowed_cmds = ReseqTrack::Tools::RunPicard->get_valid_commands;
throw(
  "Don't recognise command $command. Acceptable commands are: @allowed_cmds")
  if ( !grep { $command eq $_ } @allowed_cmds );

my @allowed_options = keys %{&ReseqTrack::Tools::RunPicard::DEFAULT_OPTIONS};
foreach my $option ( keys %options ) {
  throw(
    "Don't recognise option $option. Acceptable options are: @allowed_options")
    if ( !grep { $option eq $_ } @allowed_options );
}

throw("Must specify an output directory") if ( !$output_dir );
throw("Must specify an output type")      if ( !$type_output );
throw("Must specify an index type if create_index flag is used")
  if ( $create_index && !$type_index );

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

my $collection = $ca->fetch_by_name_and_type( $name, $type_input );
throw("Failed to find a collection for "
    . $name . " "
    . $type_input
    . " from "
    . $dbname )
  unless ($collection);

my $input_files = $collection->others;
my @input_filepaths = map { $_->{'name'} } @$input_files;

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
    $output_dir =
      create_directory_path( $run_meta_info, $directory_layout, $output_dir );
  }
}

my $picard_object = ReseqTrack::Tools::RunPicard->new(
  -input_files  => \@input_filepaths,
  -working_dir  => $output_dir,
  -job_name     => $name,
  -options      => \%options,
  -java_exe     => $java_exe,
  -jvm_options  => $jvm_options,
  -picard_dir   => $picard_dir,
  -create_index => $create_index,
  -keep_metrics => (defined $metrics_file_type),
);
my ($metrics) = $picard_object->run($command);

throw("store_stats is set, but the command has not produced any")  if ( $store_stats && !$metrics );

$db->dbc->disconnect_when_inactive(0);
if ($store) {
  my $host = get_host_object( $host_name, $db );

  my @file_paths = @{ $picard_object->output_bam_files };

  if (@file_paths) {
    foreach my $path (@file_paths) {
      throw("database already has file with name $path")
        if ( $fa->fetch_by_name($path) );
    }
    my $files =
      create_objects_from_path_list( \@file_paths, $type_output, $host );
    if ( !$disable_md5 ) {
      foreach my $file (@$files) {
        $file->md5( run_md5( $file->name ) );
      }
    }


    my $collection = ReseqTrack::Collection->new(
      -name   => $name,
      -type   => $type_output,
      -others => $files
    );
    
    $ca->store($collection);
 
    if ( $store_stats && $metrics ) {
      for my $metrics_row (@$metrics) {
        while ( my ( $key, $value ) = each %$metrics_row ) {
          if ( defined $value && defined $key ) {
            my $stat = create_statistic_for_object( $collection, $key, $value );
            $collection->statistics($stat);
          }
        }
      }
      $ca->store_statistics($collection);
    }

    if ( $metrics_file_type && @{ $picard_object->output_metrics_files } ) {
      create_metrics_records( $name, $metrics_file_type,
        $picard_object->output_metrics_files );
    }
  }

  my $fa        = $db->get_FileAdaptor;
  my $bai_paths = $picard_object->output_bai_files;
  if (@$bai_paths) {
    my $bais = create_objects_from_path_list( $bai_paths, $type_index, $host );
    foreach my $bai (@$bais) {
      if ( !$disable_md5 ) {
        $bai->md5( run_md5( $bai->name ) );
      }
      $fa->store($bai);
    }
  }
}

if ($delete_inputs) {
  my @index_files;
  my @index_filepaths;
INDEX_FILE:
  foreach my $filepath (@input_filepaths) {
    my $index_filepath = $filepath . '.bai';
    next INDEX_FILE if ( !-f $index_filepath );
    push( @index_filepaths, $index_filepath );
    my $index_file = $fa->fetch_by_name($index_filepath);
    if ($index_file) {
      push( @index_files, $index_file );
    }
  }

  my $ha = $db->get_HistoryAdaptor;
  foreach my $file ( @$input_files, @index_files ) {
    my $history = ReseqTrack::History->new(
      -other_id   => $file->dbID,
      -table_name => 'file',
      -comment    => "deleted by $0"
    );
    $ha->store($history);
  }
  foreach my $file ( @input_filepaths, @index_filepaths ) {
    delete_file($file);
  }
}

sub create_metrics_records {
  my ( $collection_name, $type, $file_names ) = @_;

  #create a File loader object for the files
  my $loader = ReseqTrack::Tools::Loader::File->new(
    -file            => $file_names,
    -do_md5          => 1,
    -hostname        => $host_name,
    -db              => $db,
    -assign_types    => 0,
    -check_types     => 0,
    -type            => $type,
    -update_existing => 1,
  );

  #Load the files into the database
  $loader->process_input();
  $loader->create_objects();
  $loader->sanity_check_objects();
  my $file_objects = $loader->load_objects();

  my $collection = ReseqTrack::Collection->new(
    -name       => $collection_name,
    -others     => $file_objects,
    -type       => $type,
    -table_name => 'file',
  );

  $db->get_CollectionAdaptor->store($collection,1);
}

=pod

=head1 NAME

reseqtrack/scripts/process/run_picard.pl

=head1 SYNOPSIS

This script runs picard to process sam / bam files.  It will do one of the following;
  
      merge: all sam / bam files in the collection are merged into a single sorted bam file
      mark_duplicates: mark or delete duplicates for each sam / bam file in a collection
      fix_mate: fix mate information for each sam / bam file in a collection
      sort: sorts each sam / bam file in a collection (does not merge)
      alignment_metrics: make a metrics file for each sam / bam file in a collection

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

  -picard_dir, path to the directory containing the picard jar files
  NB the RunPicard class can guess the location of picard_dir if the PICARD environment variable is set

  -java_exesamtools, the java executable.  NB the RunPicard class uses the default 'java'

  -jvm_args, options of java. The RunPicard class uses default values if nothing is specified.

  -host_name, default is '1000genomes.ebi.ac.uk', needed for storing output files

  -store, boolean flag, to store output files in the database.

  -disable_md5, boolean flag, files written to the database will not have an md5

  -delete_inputs, boolean flag, to delete the original input files (and index if it exists)
  and remove mark them as deleted in the database History table

  -directory_layout, specifies where the files will be located under output_dir.
      Tokens matching method names in RunMetaInfo will be substituted with that method's
      return value.

  -run_id_regex, used to get run meta info.  Default is '[ESD]RR\d{6}'
  -study_id_regex, used to get run meta info.  Default is '[ESD]RS\d{6}'

  -command, must be one of the following: 'merge', 'sort', 'mark_duplicates', 'alignment_metrics', 'fix_mate'
  Tells the RunPicard object how to process the input files

  -options, for constructing the options hash passed to the RunPicard object
  e.g. -options assume_sorted=1 -options index_ext=.bai

  -create_index, flag to create an index file for any output bam file

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/process/run_picard.pl $DB_OPTS -command merge -create_index -name SRR000001
    -type_input SAM -type_output MERGED_BAM -type_index MERGED_BAI
    -store -output_dir /path/to/base_dir/ -directory_layout population/sample_id/run_id
    -options use_threading=1 -options validation_stringency=SILENT

=cut

