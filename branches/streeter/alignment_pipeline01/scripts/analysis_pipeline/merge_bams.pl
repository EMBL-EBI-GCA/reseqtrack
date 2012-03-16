#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_object_from_path);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 delete_file);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use ReseqTrack::Tools::RunSamtools;
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
my $index_bams = 1;
my $directory_layout;

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
  'index_bams=s' => \$index_bams,
  'directory_layout=s' => \$directory_layout,
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
my $fa = $db->get_FileAdaptor;

my $collection = $ca->fetch_by_name_and_type($name, $type_input);
throw("Failed to find a collection for ".$name." ".$type_input." from ".$dbname) 
    unless($collection);

my $input_files = $collection->others;
my @input_filepaths = map {$_->{'name'}} @$input_files;

throw("no output_dir") if (!$output_dir);
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
                    -output_to_working_dir   => 1,
                    -flag_merge              => 1,
                    -flag_index              => $index_bams,
                    );
$samtools_object->run;

if($store){
  my $host = get_host_object($host_name, $db);

  my $output_bams = $samtools_object->output_bam_files;
  throw("found " . @$output_bams . " bams\n") if (@$output_bams != 1);
  my $bam_path = $output_bams->[0];
  throw("database already has file with name $bam_path")
      if ($fa->fetch_by_name($bam_path));
  my $bam = create_object_from_path($bam_path, $type_output, $host);
  $bam->md5( run_md5($bam->name) );
  my $collection = ReseqTrack::Collection->new(
      -name => $name, -type => $type_output,
      -others => $bam);
  $ca->store($collection);

  if ($index_bams) {
    my $fa = $db->get_FileAdaptor;
    my $bai_paths = $samtools_object->output_bai_files;
    throw("found " . @$bai_paths . " bais\n") if (@$bai_paths != 1);
    my $bai = create_object_from_path($bai_paths->[0], $type_index, $host);
    $bai->md5( run_md5($bai->name) );
    $fa->store($bai);
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
