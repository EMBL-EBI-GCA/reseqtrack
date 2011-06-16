#!/sw/arch/bin/perl -w

use strict;
use Getopt::Long;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Collection;
use ReseqTrack::Host;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::ERAUtils;
use File::Basename;
use File::Path;

$| = 1;

my $run_id;
my $ftp_server = 'ftp://ftp.sra.ebi.ac.uk';
my $program = 'wget';
my $options = '-t 1 --ignore-length';
my $root_dir;
my $host_name;
my $file_type;
my $clobber;
my $remote;
my $dbhost;
my $dbname;
my $dbuser;
my $dbport;
my $dbpass;
my $output_dir;
my $load = 0;
my $era_dbuser;
my $era_dbpass;
my $help = 0;

&GetOptions( 
  'run_id=s' => \$run_id,
  'program=s' => \$program,
  'options=s' => \$options,
  'output_dir=s' => \$output_dir,
  'host_name=s' => \$host_name,
  'type|file_type=s' => \$file_type,
  'clobber!' => \$clobber,
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'remote!' => \$remote,
  'load!' => \$load,
  'era_dbuser=s' =>\$era_dbuser,
  'era_dbpass=s' => \$era_dbpass,
  'help!' => \$help,
    );

if($help){
  useage();
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $era_db = get_erapro_conn($era_dbuser, $era_dbpass);

my $rmi_a = $db->get_RunMetaInfoAdaptor;
my $meta_info = $rmi_a->fetch_by_run_id($run_id);
throw("Failed to find a run meta info object for ".$run_id." from ".$dbname)
    unless($meta_info);
if($meta_info->status ne 'public'){
  exit(0);
}
my $full_output_dir = $output_dir."/".$meta_info->sample_name."/archive_sequence";
$full_output_dir =~ s/\/\//\//;
unless(-d $full_output_dir){
  mkpath($full_output_dir);
}

my $era_rmia = $era_db->get_ERARunMetaInfoAdaptor;

unless($era_rmia->is_fastq_available($run_id)){
  print STDERR "There are no fastq files available for ".$run_id."\n";
  exit(20);
}
$db->dbc->disconnect_when_inactive(1);
my $hash = get_era_fastq($run_id, $full_output_dir, $ftp_server, $clobber, $program, $options, $era_db);


my $ha = $db->get_HostAdaptor;
my $host = $ha->fetch_by_name($host_name);
if(!$host){
  $host = ReseqTrack::Host->new
      (
       -name => $host_name,
       -remote => $remote
      );
}

my @objects;

throw("Have failed to fetch any fastq files for ".$run_id." in ".$full_output_dir)
    unless(keys(%$hash));


foreach my $path(keys(%$hash)){
  print $path."\n";
  my $size = -s $path;
  if($size <= 20){
    throw($path." appears to be empty");
  }
  my $files = create_objects_from_path_list([$path], $file_type, $host);
  my $file = $files->[0];
  my $md5 = $hash->{$file->name};
  if(!$md5){
    warning("Don't have an md5 for ".$file->filename);
  }
  $file->md5($md5);
  push(@objects, $file);
}
$db->dbc->disconnect_when_inactive(0);
if($load){
  my $fa = $db->get_FileAdaptor;
  foreach my $object(@objects){
    my $objects = $fa->fetch_by_filename($object->filename);
    if($objects && @$objects >= 1){
      my $stored_object;
      if(@$objects == 1){
        $stored_object = $objects->[0];
      }else{
        throw("Not sure what to do, seem to have ".@$objects." associated with ".
              $object->filename);
      }
      my $history = create_history($object, $stored_object);
      $object->dbID($stored_object->dbID);
      $object->history($history);
    }
    $fa->store($object, 1);
  }
}

my $collection =  ReseqTrack::Collection->new(
  -name => $run_id,
  -others => \@objects,
  -type => $file_type,
  -table_name => 'file',
    );

my $ca = $db->get_CollectionAdaptor;
$ca->store($collection) if($load);


my $new_collection = $ca->fetch_by_name_and_type($run_id, $file_type);
my $others = $new_collection->others;
my @paths;
foreach my $other(@$others){
  push(@paths, $other->name);
  throw($other->name." doesn't exist") unless(-e $other->name);
}

sub useage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/process/get_fastq.pl

=head1 SYNOPSIS

This script fetches fastq files from the ERA ftp servers.

=head1 OPTIONS

Database options

This is the database the objects will be written to

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the standard
port for mysql-g1kdcc.ebi.ac.uk

These are for the ERA database

-era_dbuser, the name of the era database user
-era_dbpass, the password for the era database

Input options

-run_id, this is the run id for the run you wish to fetch the fastq files for
-program, this is the program to use to run the fetch, by default this is wget
-options, this is any commandline options for the above program
-output_dir, this is the root path where the files should be written. Meta information associated with the run id will be used to complete the path in the style root_path/sample_id/archive_fastq
-type/file_type, this is the file/collection type to use when loading the data
-host_name, this is the name of the host object to associated with the files
-clobber, this specifies whether existing files should be copied over or not
-remote, this specifies if the give host is remote
-load, this species if the created file/collection objects should be loaded into
the database

=head1 Examples

perl ReseqTrack/scripts/process/get_fastq.pl -output scratch/staging_area/sequence_staging/ -clobber -dbhost host -dbuser user -dbpass **** -dbport 4197 -dbname staging_db -host 1000genomes.ebi.ac.uk -file_type ARCHIVE_FASTQ -load -run_id SRR031624 -era_dbuser username -era_dbpass *****

=cut
