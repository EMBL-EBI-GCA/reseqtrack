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
my $delete_input;
my $realign;
my $recalibrate;

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
  'delete_input!' => \$delete_input,
  'realign!' => \$realign,
  'recalibrate!' => \$recalibrate,
    );

throw("must specify one of -realign or -recalibrate, not both or neither")
    if (($realign && $recalibrate) || (!$realign && !$recalibrate));
my $gatk_module = load_gatk_module($realign ? 'IndelRealigner' : 'QualityScoreRecalibrator');

my @allowed_options = keys %{eval "$gatk_module".'::DEFAULT_OPTIONS'};
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
                  -samtools                => $samtools,
                  -job_name                => $name,
                  -java_exe                => $java_exe,
                  -jvm_args                => $jvm_args,
                  -gatk_path               => $gatk_path,
                  -options                 => \%options,
                  -known_sites_files      => \@known_sites_vcf,
                  );
$gatk_object->run;

if($store){
  my $host = get_host_object($host_name, $db);

  my $bam = create_object_from_path($gatk_object->output_bam, $type_output, $host);
  $bam->md5( run_md5($bam->name) );
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
