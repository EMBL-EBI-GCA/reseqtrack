#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 );
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use Getopt::Long;

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
my $module_path = 'ReseqTrack::Tools::RunAlignment';
my $module_name;
my %module_constructor_args;
my $program_file;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $update;
my $type_output = 'SAM';
my $paired_length;
my $first_read;
my $last_read;
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
  'module_path=s' => \$module_path,
  'module_name=s' => \$module_name,
  'module_constructor_arg=s' => \%module_constructor_args,
  'program_file=s' => \$program_file,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'update!' => \$update,
  'paired_length=s' => \$paired_length,
  'first_read=s' => \$first_read,
  'last_read=s' => \$last_read,
  'RG_field=s' => \%read_group_fields,
  'directory_layout=s' => \$directory_layout,
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
my $ca = $db->get_CollectionAdaptor;

my $collection = $ca->fetch_by_name_and_type($name, $type_input);

throw("Failed to find a collection for ".$name." ".$type_input." from ".$dbname) 
    unless($collection);

my $input_files = $collection->others;
my @input_file_paths = map {$_->{'name'}} @$input_files;

if ($name =~ /[ESD]RR\d{6}/) {
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

if (! defined $first_read && ! defined $last_read) {
  if ($name =~ /_(?:FRAG|MATE)(\d+),(\d+)/i) {
    $first_read = $1;
    $last_read = $2;
  }
}

my $alignment_module = $module_path."::".$module_name;
my $constructor_hash;
while (my ($key, $value) = each %module_constructor_args) {
  $constructor_hash->{'-'.$key} = $value;
}
$constructor_hash->{-input_files} = \@input_file_paths;
$constructor_hash->{-program} = $program_file;
$constructor_hash->{-working_dir} = $output_dir;
$constructor_hash->{-reference} = $reference;
$constructor_hash->{-job_name} = $name;
$constructor_hash->{-paired_length} = $paired_length;
$constructor_hash->{-first_read} = $first_read;
$constructor_hash->{-last_read} = $last_read;
$constructor_hash->{-read_group_fields} = \%read_group_fields;

my $run_alignment = setup_alignment_module($alignment_module, $constructor_hash);

$run_alignment->run;


if($store){
  my $host = get_host_object($host_name, $db);
  my $fa = $db->get_FileAdaptor;
  my $ca = $db->get_CollectionAdaptor;

  my $sam_filepaths = $run_alignment->output_files;
  my $sam_files = create_objects_from_path_list($sam_filepaths, $type_output, $host);
  foreach my $file(@$sam_files){
    my $md5 = run_md5($file->name);
    $file->md5($md5);
    $fa->store($file, $update);
  }
  my $collection = ReseqTrack::Collection->new(
      -name => $name, -type => $type_output,
      -others => $sam_files);
  $ca->store($collection);
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
