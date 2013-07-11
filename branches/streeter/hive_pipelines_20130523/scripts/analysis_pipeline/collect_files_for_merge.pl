#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::FileSystemUtils qw(delete_file check_directory_exists);
use ReseqTrack::Tools::FileUtils qw(create_object_from_path);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use Getopt::Long;
use File::Copy qw(move);

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $host_name = '1000genomes.ebi.ac.uk';
my $type_key;
my @type_input;
my $type_collection;
my $type_merged;
my $level;
my $output_dir;
my $directory_layout;
my $help;
my $verbose;
my $run_id_regex = '[ESD]RR\d{6}';
my $sample_id_regex = '[ESD]RS\d{6}';

&GetOptions( 
    'dbhost=s'      => \$dbhost,
    'dbname=s'      => \$dbname,
    'dbuser=s'      => \$dbuser,
    'dbpass=s'      => \$dbpass,
    'dbport=s'      => \$dbport,
    'host_name=s' => \$host_name,
    'level=s'      => \$level,
    'type_key=s'       => \$type_key,
    'type_input=s'       => \@type_input,
    'type_collection=s'       => \$type_collection,
    'type_merged=s'       => \$type_merged,
    'output_dir=s' => \$output_dir,
    'directory_layout=s' => \$directory_layout,
    'help!'         => \$help,
    'verbose' => \$verbose,
    'run_id_regex=s' => \$run_id_regex,
    'sample_id_regex=s' => \$sample_id_regex,
   );


if ($help) {
    exec('perldoc', $0);
    exit(0);
}

throw("Must specify an output directory") if (!$output_dir);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
                         -host => $dbhost,
                         -user => $dbuser,
                         -port => $dbport,
                         -dbname => $dbname,
                         -pass => $dbpass,
                        );

my $ca = $db->get_CollectionAdaptor;
my $rmia = $db->get_RunMetaInfoAdaptor;

my %inputs;
foreach my $type_input (@type_input) {
  foreach my $collection (@{$ca->fetch_by_type($type_input)}) {
    $inputs{$collection->name} = $collection->others;
    print "Found input collection: ".$collection->name.$/ if $verbose;
  }
}

my %existing_output_collections;
foreach my $collection (@{$ca->fetch_by_type($type_collection)}) {
  $existing_output_collections{$collection->name} = 1;
  print "Found output collection: ".$collection->name.$/ if $verbose;
}
my %existing_merged;
foreach my $collection (@{$ca->fetch_by_type($type_merged)}) {
  $existing_merged{$collection->name} = 1;
  print "Found merged collection: ".$collection->name.$/ if $verbose;
}

my %files_to_move;

if ($level eq 'RUN') {
  my $key_collections = $ca->fetch_by_type($type_key);

  COLLECTION:
  foreach my $key_collection (@$key_collections) {
    next COLLECTION if ($existing_output_collections{$key_collection->name});
    next COLLECTION if ($existing_merged{$key_collection->name});
    my @others;
    foreach my $name (map {$_->name} @{$key_collection->others}) {
      next COLLECTION if (!exists $inputs{$name});
      push(@others, @{$inputs{$name}});
    }
    if (@others == 1) {
      $files_to_move{$key_collection->name} = $others[0];
      print "Found file to move for colleciton: ".$key_collection->name.$/ if $verbose;
    }
    else {
      my $output_collection = ReseqTrack::Collection->new(
              -name => $key_collection->name, -type => $type_collection,
              -others => \@others);
      $ca->store($output_collection);
      print "Created merge collection: ".$key_collection->name.$/ if $verbose;
    }
  }
}


if ($level eq 'LIBRARY') {
  my %library_runs;
  foreach my $rmi (@{$rmia->fetch_all}) {
    push(@{$library_runs{$rmi->sample_id}{$rmi->library_name}}, $rmi->run_id);
  }
  foreach my $sample_id (keys %library_runs) {
    LIBRARY:
    while (my ($library_name, $run_id_list) = each %{$library_runs{$sample_id}}) {
      my $output_name = $sample_id.'_'.$library_name;
      next LIBRARY if ($existing_output_collections{$output_name});
      next LIBRARY if ($existing_merged{$output_name});
      my @others;
      foreach my $run_id (@$run_id_list) {
        next LIBRARY if (!exists $inputs{$run_id});
        push(@others, @{$inputs{$run_id}});
      }
      if (@others == 1) {
        $files_to_move{$output_name} = $others[0];
        print "Found file to move for library: $library_name $/" if $verbose;
      }
      else {
        my $output_collection = ReseqTrack::Collection->new(
                -name => $output_name, -type => $type_collection,
                -others => \@others);
        $ca->store($output_collection);
        print "Created merge library: $library_name $/" if $verbose;
      }
    }
  }
}

if ($level eq 'SAMPLE') {
  my %sample_libraries;
  foreach my $rmi (@{$rmia->fetch_all}) {
    $sample_libraries{$rmi->sample_id}{$rmi->library_name} = 1;
  }
  SAMPLE:
  foreach my $sample_id (keys %sample_libraries) {
    next SAMPLE if ($existing_output_collections{$sample_id});
    next SAMPLE if ($existing_merged{$sample_id});
    my @others;
    foreach my $library_name (keys %{$sample_libraries{$sample_id}}) {
      my $inputs_key = $sample_id.'_'.$library_name;
      next SAMPLE if (!exists $inputs{$inputs_key});
      push(@others, @{$inputs{$inputs_key}});
    }
    if (@others == 1) {
      $files_to_move{$sample_id} = $others[0];
      print "Found file to move for sample: $sample_id $/" if $verbose;
    }
    else {
      my $output_collection = ReseqTrack::Collection->new(
              -name => $sample_id, -type => $type_collection,
              -others => \@others);
      $ca->store($output_collection);
      print "Created merge sample: $sample_id $/" if $verbose;
    }
  }
}


my $ha = $db->get_HistoryAdaptor;
my $fa = $db->get_FileAdaptor;
my $host = get_host_object($host_name, $db);
while (my ($name, $file) = each %files_to_move) {

  my $file_output_dir = $output_dir;
  if ($directory_layout) {
    my $run_meta_info;
    if ($name =~ /$run_id_regex/) {
      $run_meta_info = $rmia->fetch_by_run_id($&);
    }
    elsif ($name =~ /$sample_id_regex/) {
      my $rmi_list = $rmia->fetch_by_sample_id($&);
      $run_meta_info = $rmi_list->[0] if (@$rmi_list);
    }

    if ($run_meta_info) {
      $file_output_dir = create_directory_path($run_meta_info, $directory_layout, $output_dir);
    }
  }
  check_directory_exists($file_output_dir);

  my $old_file_path = $file->name;
  my $new_file_path = "$file_output_dir/$name.bam";

  my $history = ReseqTrack::History->new(
      -other_id => $file->dbID, -table_name => 'file',
      -comment => "moved by $0 to $new_file_path");
  move ($old_file_path, $new_file_path)
    or throw ("could not move $old_file_path to $new_file_path $!");
  $ha->store($history);

  my $old_index_path = $old_file_path . '.bai';
  if (-f $old_index_path) {
    delete_file($old_index_path);
    my $old_index = $fa->fetch_by_name($old_index_path);
    if ($old_index) {
      my $index_history = ReseqTrack::History->new(
        -other_id => $old_index->dbID, -table_name => 'file',
        -comment => "deleted by $0");
      $ha->store($history);
    }
  }

  my $merged_file = create_object_from_path($new_file_path, $type_merged, $host);
  $merged_file->md5( $file->md5 );
  my $output_collection = ReseqTrack::Collection->new(
          -name => $name, -type => $type_merged,
          -others => [$merged_file]);
  $ca->store($output_collection);

}

=pod

=head1 NAME

reseqtrack/scripts/event/collect_files_for_merge.pl

=head1 SYNOPSIS

This script finds files that are ready to be merged for the next level of processing
(i.e. run level, library level or sample level)
Files must stored in collections, whose names are used to determine how they will be
grouped together in the output collections.
If a collection contains only one file to be 'merged' then it will be moved to a new
location and given a new type (i.e. it mimics the effects of a merge event)

For merging to sample level, the input collections must be named SampleID_LibraryName.
The output collection will be named by the sample_id.

For merging to library level, the input collections must be named by the run_id.
The output collection will be named SampleID_LibraryName.

For merging to run level, a 'key_collection' will be used to determine the names of the
input collections. The name of key_collection must be run_id, and the names of
key_collection->others will be the identifiers. Usually the output collection from
make_split_fastq_collections.pl should be used as the key_collection.
The output collection will be named by the run_id

=head1 OPTIONS

  database options:

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on

  other options:

  -level, must be either 'RUN', 'LIBRARY' or 'SAMPLE'

  -type_input, type of the collections that will be searched for files. Can be specified
  multiple times

  -type_collection, type of the output collection that will be written for files that
  need to be merged

  -type_key, type of the 'key_collection' as described above.  Only needed when merging
  to run level

  -type_merged, type of the output collection that will be written for files that do not
  need to be merged

  -output_dir, base directory to hold files that do not need to be merged

  -directory_layout, specifies where the files will be located under output_dir.
      Tokens matching method names in RunMetaInfo will be substituted with that method's
      return value.

  -run_id_regex, used to get run meta info.  Default is '[ESD]RR\d{6}'
  -study_id_regex, used to get run meta info.  Default is '[ESD]RS\d{6}'

  -help, flag to print this help and exit

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/event/collect_files_for_merge.pl $DB_OPTS -level SAMPLE
    -type_input LIBRARY_BAM -type_collection LIBRARY_BAM_SET
    -type_merged SAMPLE_BAM
    -output_dir /path/to/base/dir -directory_layout population/sample_id

=cut

