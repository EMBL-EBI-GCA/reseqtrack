#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 );
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
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
my %options;

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
  'option=s' => \%options,
    );

my $alignment_module = load_alignment_module($module_path, $module_name);
my $allowed_options = get_allowed_options($alignment_module);
foreach my $option (keys %options) {
  throw("Don't recognise option $option. Acceptable options are: ".join(' ', @$allowed_options))
    if (! grep {$option eq $_ } @$allowed_options);
}

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
  if ($name =~ /_(?:FRAG|MATE)(\d+)-(\d+)/i) {
    $first_read = $1;
    $last_read = $2;
  }
}

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
$constructor_hash->{-options} = \%options;

my $run_alignment = $alignment_module->new(%$constructor_hash);

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

sub load_alignment_module {
  my ($module_path, $module_name) = @_;
  my $full_module_path = join('::', $module_path, $module_name);
  my $file = "$full_module_path.pm";
  $file =~ s{::}{/}g;
  eval {
    require "$file";
  };
  if ($@) {
    throw("cannot load $file: $@")
  }
  return $full_module_path;
}

sub get_allowed_options {
  my $alignment_module = shift;
  my $default_options_sub = '&'.$alignment_module . '::DEFAULT_OPTIONS';
  return [] if (! eval "exists $default_options_sub");
  my $default_options = eval "$default_options_sub";
  my @allowed_options = keys %$default_options;
  return \@allowed_options;
}

=pod

=head1 NAME

reseqtrack/scripts/process/run_alignment.pl

=head1 SYNOPSIS

This script runs an alignment using any child class of ReseqTrack::Tools::RunAlignment

=head1 OPTIONS

  database options:

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on

  other options:

  -name, name of the collection of fastq files
  If name is a run_id (or contains a run_id), the run_meta_info table will be used to get some info

  -type_input, type of the collection of fastq files

  -output_dir, base directory to hold files that do not need to be merged

  -reference, path to the reference fasta file for the aligner

  -module_path, the default is 'ReseqTrack::Tools::RunAlignment'

  -module_name, should be e.g. 'Stampy', 'BWA', 'BFAST', 'Smalt'

  -module_constructor_arg, command line format is -module_constructor_arg arg_name=arg_value.
  Can be specified multiple times on the command line.
  These constructor args are passed directly to the RunAlignment object.

  -program_file, executable of the alignment program
  NB some of the RunAlignment child classes can guess at its location if it is not specified

  -store, boolean flag, to store output files in the database.

  -host_name, default is '1000genomes.ebi.ac.uk', needed for storing output files

  -type_output, collection type and file type when storing output files in the database

  -paired_length, integer, insert size passed to the RunAlignment object.
  Will be read from the run_meta_info table if not specified on the command line

  -first_read, integer, start aligning from this read in the fastq file.
  If not specified, it will be inferred from 'name' if name is something like ERR000001_FRAG1000-2000

  -last_read, integer, stop aligning at this read in the fastq file.
  If not specified, it will be inferred from 'name' if name is something like ERR000001_FRAG1000-2000

  -RG_field, e.g. -RG_field ID=ERR000001 -RG_field LB=my_lib -RG_field SM=my_sample
  If no RG fields are specified they will be taken from the run_meta_info table

  -directory_layout, specifies where the files will be located under output_dir.
      Tokens matching method names in RunMetaInfo will be substituted with that method's
      return value.

  -option, for constructing the options hash passed to the RunAlignment object
  e.g. -option threads=4 -option

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/process/run_alignment.pl $DB_OPTS -name SRR000001
    -type_input FASTQ -type_output SAM -reference /path/to/ref
    -module_name BWA -option threads=4 -program_name /path/to/bwa
    -store -output_dir /path/to/base_dir/ -directory_layout population/sample_id/run_id

=cut

