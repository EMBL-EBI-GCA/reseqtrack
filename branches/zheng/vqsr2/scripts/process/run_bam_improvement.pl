#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_object_from_path);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 delete_file);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
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
my $output_dir;
my $java_exe;
my $jvm_args;
my $gatk_path;
my $reference;
my %options;
my @known_sites_vcf;
my $samtools;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $disable_md5;
my $delete_input;
my $realign;
my $recalibrate;
my $directory_layout;
my $intervals_file;

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
  'java_exe=s' => \$java_exe,
  'jvm_args=s' => \$jvm_args,
  'gatk_path=s' => \$gatk_path,
  'reference=s' => \$reference,
  'options=s' => \%options,
  'known_sites_vcf=s' => \@known_sites_vcf,
  'samtools=s' => \$samtools,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'disable_md5!' => \$disable_md5,
  'delete_input!' => \$delete_input,
  'realign!' => \$realign,
  'recalibrate!' => \$recalibrate,
  'directory_layout=s' => \$directory_layout,
  'intervals_file=s' => \$intervals_file,
    );

throw("must specify one of -realign or -recalibrate, not both or neither")
    if (($realign && $recalibrate) || (!$realign && !$recalibrate));
my $gatk_module = load_gatk_module($realign ? 'IndelRealigner' : 'QualityScoreRecalibrator');

my @allowed_options = keys %{eval '&'."$gatk_module".'::DEFAULT_OPTIONS'};
foreach my $option (keys %options) {
  throw("Don't recognise option $option. Acceptable options are: @allowed_options")
    if (! grep {$option eq $_ } @allowed_options);
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

my $others = $collection->others;
throw("Expecting one file; have ".@$others."\n")
    if (@$others != 1);
my $input_file = $others->[0];

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

my $gatk_object = $gatk_module->new(
                  -input_files             => $input_file->name,
                  -working_dir             => $output_dir,
                  -reference               => $reference,
                  -program                 => $samtools,
                  -job_name                => $name,
                  -java_exe                => $java_exe,
                  -jvm_args                => $jvm_args,
                  -gatk_path               => $gatk_path,
                  -options                 => \%options,
                  -known_sites_files      => \@known_sites_vcf,
                  );
if ($realign) {
  $gatk_object->intervals_file($intervals_file);
}
$gatk_object->run;

if($store){
  my $host = get_host_object($host_name, $db);

  my $bam = create_object_from_path($gatk_object->output_bam, $type_output, $host);
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


sub load_gatk_module {
  my $module_name = shift;
  my $file = "ReseqTrack/Tools/GATKTools/$module_name.pm";
  eval {
    require "$file";
  };
  if ($@) {
    throw("cannot load $file: $@")
  }
  my $module = "ReseqTrack::Tools::GATKTools::$module_name";
  return $module;
}

=pod

=head1 NAME

reseqtrack/scripts/process/run_bam_improvement.pl

=head1 SYNOPSIS

This script is used to run either indel realignment or quality score recalibration on a bam file

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

  -realign, boolean flag, run the Indel Realigner

  -recalibrate, boolean flag, run the Quality Score Recalibrator

  -name, name of a collection containing the bam file
  If name is a run_id (or contains a run_id), the run_meta_info table will be used to get some info

  -type_input, type of the collection containing the input bam file

  -type_output, collection type and file type when storing output bam file in the database

  -output_dir, base directory to hold files that do not need to be merged

  -java_exe, path to the java executable. The GATKTools class can guess at its location if it is not specified.

  -jvm_args, options of java. The GATKTools class uses default values if nothing is specified.

  -gatk_path, path to a directory containing the gatk jar files.
  Alternatively, the GATKTools class will use the $GATK environment variable if set.

  -reference, path to the reference fasta file

  -options, for constructing the options hash passed to the GATKTools object
  e.g. -options threads=4

  -known_sites_vcf, path to a vcf file.  Can be specified multiple times.
  For the Indel Realigner, this contains known indels
  For the Quality Score Recalibrator, this contains known polymorphic sites

  -intervals_file, path to an intervals file, only needed for the Indel Realigner with the knowns_only option set

  -samtools, path to samtools executable. This is needed for creating index files
  NB the Samtools class can guess the location of samtools if not specified.

  -store, boolean flag, to store output bam file in the database.

  -disable_md5, boolean flag, files written to the database will not have an md5

  -host_name, default is '1000genomes.ebi.ac.uk', needed for storing output bam file

  -delete_input, boolean flag, to delete the original bam file (and index if it exists)
  and remove mark them as deleted in the database History table

  -directory_layout, specifies where the files will be located under output_dir.
      Tokens matching method names in RunMetaInfo will be substituted with that method's
      return value.

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/process/run_bam_improvement.pl $DB_OPTS -realign -name SRR000001
    -type_input BAM -type_output REALIGNED_BAM -output_dir /path/to/base_dir/
    -reference /path/to/ref -options threads=4 -known_sites_vcf /path/to/vcf
    -store -delete_input -directory_layout population/sample_id/run_id

=cut

