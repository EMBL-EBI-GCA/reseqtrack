#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 delete_file);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use File::Basename qw(fileparse);
use Getopt::Long;
use List::Util qw(first);

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $name;
my $type_input;
my $output_dir;
my $reference;
my $ref_samtools_index;
my $module_path = 'ReseqTrack::Tools::RunAlignment';
my $module_name;
my $module_constructor_args;
my $program_file;
my $samtools;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $update;
my $delete_inputs;
my $type_output = 'BAM';
my $type_index = 'BAI';
my $merge_bams;
my $sort_bams;
my $index_bams;
my $convert_sam_to_bam;
my $paired_length;
my %read_group_fields;
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
  'reference=s' => \$reference,
  'ref_samtools_index=s' => \$ref_samtools_index,
  'module_path=s' => \$module_path,
  'module_name=s' => \$module_name,
  'module_constructor_args=s' => \$module_constructor_args,
  'program_file=s' => \$program_file,
  'samtools=s' => \$samtools,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'delete_inputs!' => \$delete_inputs,
  'update!' => \$update,
  'merge_bams!' => \$merge_bams,
  'index_bams!' => \$index_bams,
  'sort_bams!' => \$sort_bams,
  'convert_sam_to_bam!' => \$convert_sam_to_bam,
  'paired_length=s' => \$paired_length,
  'RG_field=s' => \%read_group_fields,
  'directory_layout=s' => \$directory_layout,
    );

$sort_bams ||= $merge_bams;
$convert_sam_to_bam ||= $sort_bams || $index_bams;

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

my $input_files = $collection->others;
my @input_file_paths = map {$_->{'name'}} @$input_files;

if (!$output_dir) {
  $output_dir = (fileparse($input_file_paths[0]))[1];
  my $file_type = $input_files->[0]->type;
  $output_dir =~ s/$file_type\/*$/$type_output/;
}

if (!$paired_length) {
  foreach my $input_file (@$input_files) {
    my $statistic = first {$_->attribute_name eq 'paired_length'} @{$input_file->statistics};
    if ($statistic) {
      $paired_length = $statistic->attribute_value;
    }
    last if ($paired_length);
  }
}

if ($name =~ /[ESD]RR\d{6}/ &&
      ($directory_layout || !$paired_length || !%read_group_fields)) {
  my $rmia = $db->get_RunMetaInfoAdaptor;
  my $run_meta_info = $rmia->fetch_by_run_id($&);
  if ($run_meta_info && $directory_layout) {
    $output_dir = create_directory_path($run_meta_info, $directory_layout, $output_dir);
  }
  if ($run_meta_info && !$paired_length) {
    $paired_length = $run_meta_info->paired_length;
  }
  if (! %read_group_fields) {
    $read_group_fields{'ID'} = $run_meta_info->run_id;
    $read_group_fields{'CN'} = $run_meta_info->center_name;
    $read_group_fields{'LB'} = $run_meta_info->library_name;
    $read_group_fields{'PI'} = $run_meta_info->paired_length;
    $read_group_fields{'SM'} = $run_meta_info->sample_name;
    my $platform = $run_meta_info->instrument_platform;
    $platform =~ s/ABI_SOLID/SOLID/;
    $read_group_fields{'PL'} = $platform;
  }


}

my $alignment_module = $module_path."::".$module_name;
my $constructor_hash = parameters_hash($module_constructor_args);
$constructor_hash->{-input_files} = \@input_file_paths;
$constructor_hash->{-program} = $program_file;
$constructor_hash->{-working_dir} = $output_dir;
$constructor_hash->{-reference} = $reference;
$constructor_hash->{-ref_samtools_index} = $ref_samtools_index;
$constructor_hash->{-samtools} = $samtools;
$constructor_hash->{-job_name} = $name;
$constructor_hash->{-paired_length} = $paired_length;
$constructor_hash->{-read_group_fields} = \%read_group_fields;
$constructor_hash->{-merge_bams} = $merge_bams;
$constructor_hash->{-sort_bams} = $sort_bams;
$constructor_hash->{-index_bams} = $index_bams;
$constructor_hash->{-convert_sam_to_bam} = $convert_sam_to_bam;

my $run_alignment = setup_alignment_module($alignment_module, $constructor_hash);

$run_alignment->run;


if($store){
  my $host = get_host_object($host_name, $db);
  my $fa = $db->get_FileAdaptor;
  my $ca = $db->get_CollectionAdaptor;

  my $bam_filepaths = $run_alignment->output_bam_files;
  my $bam_files = create_objects_from_path_list($bam_filepaths, $type_output, $host);
  foreach my $file(@$bam_files){
    my $md5 = run_md5($file->name);
    $file->md5($md5);
    $fa->store($file, $update);
  }
  my $collection = ReseqTrack::Collection->new(
      -name => $name, -type => $type_output,
      -others => $bam_files);
  $ca->store($collection);

  my $bai_filepaths = $run_alignment->output_bai_files;
  if ($index_bams && @$bai_filepaths) {
    my $bai_files = create_objects_from_path_list($bai_filepaths, $type_index, $host);
    foreach my $file(@$bai_files){
      my $md5 = run_md5($file->name);
      $file->md5($md5);
      $fa->store($file, $update);
    }
  }
}

if($delete_inputs){
  my $ha = $db->get_HistoryAdaptor;
  foreach my $file (@$input_files) {
    my $history = ReseqTrack::History->new(
        -other_id => $file->dbID, -table_name => 'file',
        -comment => "deleted by $0");
    $ha->store($history);
  }
  foreach my $file (@$input_files) {
    delete_file($file->name);
  }
}

sub parameters_hash{
  my ($string) = @_;

  my %parameters_hash;

  if ($string) {
    if($string =~  /,/ || $string =~ /=>/){
      my @pairs = split (/,/, $string);
      foreach my $pair(@pairs){
        my ($key, $value) = split (/=>/, $pair);
        if ($key && ($value || $value == 0)) {
          $key   =~ s/^\s+//g;
          $key   =~ s/\s+$//g;
          $value =~ s/^\s+//g;
          $value =~ s/\s+$//g;
          $parameters_hash{$key} = $value;
        } else {
          $parameters_hash{$key} = 1;
        }
      }
    }else{
      $parameters_hash{'-options'} = $string;
    }
  }
  return \%parameters_hash;
}

sub setup_alignment_module{
  my ($alignment_module, $args) = @_;
  my $file = "$alignment_module.pm";
  $file =~ s{::}{/}g;
  eval {
    require "$file";
  };
  if($@){
    throw("ReseqTrack::Tools::EventPipeline::setup_batch_submission_system ".
          "Can't find $file [$@]");
  }
  my %constructor_args = %$args;
  my $object = $alignment_module->new(%constructor_args);
  return $object;
}
