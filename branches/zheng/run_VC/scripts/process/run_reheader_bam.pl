#!/usr/bin/env perl

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_object_from_path);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 delete_file);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);
use File::Basename qw(fileparse);
use ReseqTrack::Tools::ReheaderBam;
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
my $header_lines_file;
my $get_fastq_files;
my $type_fastq;
my $samtools;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $disable_md5;
my $delete_input;
my $directory_layout;
my $SQ_assembly;
my $SQ_species;
my $SQ_uri;
my $replace_PG;
my $replace_CO;
my $reuse_old_header;
my $root_trim;
my $ftp_root;
my $trim_paths;
my $run_id_regex = '[ESD]RR\d{6}';
my $sample_id_regex = '[ESD]RS\d{6}';

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
  'header_lines_file=s' => \$header_lines_file,
  'get_fastq_files!' => \$get_fastq_files,
  'type_fastq=s' => \$type_fastq,
  'samtools=s' => \$samtools,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'disable_md5!' => \$disable_md5,
  'delete_input!' => \$delete_input,
  'directory_layout=s' => \$directory_layout,
  'SQ_assembly=s' => \$SQ_assembly,
  'SQ_species=s' => \$SQ_species,
  'SQ_uri=s' => \$SQ_uri,
  'replace_PG!' => \$replace_PG,
  'replace_CO!' => \$replace_CO,
  'reuse_old_header!' => \$reuse_old_header,
  'root_trim=s' => \$root_trim,
  'ftp_root=s' => \$ftp_root,
  'trim_paths!' => \$trim_paths,
  'run_id_regex=s' => \$run_id_regex,
  'sample_id_regex=s' => \$sample_id_regex,
    );

throw("Must specify an output directory") if (!$output_dir);
throw("Must specify a fastq type") if ($get_fastq_files && !$type_fastq);

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



my @extra_header_lines;

if ($directory_layout || $get_fastq_files) {
  my $rmia = $db->get_RunMetaInfoAdaptor;
  my @run_meta_infos;
  if ($name =~ /$run_id_regex/) {
    my $run_meta_info = $rmia->fetch_by_run_id($&);
    push(@run_meta_infos, $run_meta_info) if $run_meta_info;
  }
  elsif ($name =~ /$sample_id_regex/) {
    my $sample_id = $&;
    my $sample_rmis = $rmia->fetch_by_sample_id($sample_id);
    if ($name =~ /${sample_id}_(\S+)/) {
      my $library_name = $1;
      @run_meta_infos = grep {$_->library_name eq $library_name} @$sample_rmis;
    }
    else {
      @run_meta_infos = @$sample_rmis;
    }
  }
  throw("did not find any run_meta_info") if (! @run_meta_infos);

  if ($directory_layout) { 
    $output_dir = create_directory_path($run_meta_infos[0], $directory_layout, $output_dir);
  }

  if ($get_fastq_files) {
    my @regexs = (qr/$run_id_regex\S*_1\.(\w+\.)*fastq(\.gz)?$/i,
                  qr/$run_id_regex\S*_2\.(\w+\.)*fastq(\.gz)?$/i,
                  qr/$run_id_regex\S*\.fastq(\.gz)?$/i);
    foreach my $rmi (@run_meta_infos) {
      my $run_id = $rmi->run_id;
      my $fastq_collection = $ca->fetch_by_name_and_type($run_id, $type_fastq);
      throw("Did not find a collection with name $run_id and type $type_fastq") if (!$fastq_collection);

      my ($mate1, $mate2, $frag) = assign_files($fastq_collection->others, \@regexs);
      my $header_fastq_type = '$' . lc($type_fastq) . '_file';
      foreach my $fastq (grep {$_} ($mate1, $mate2, $frag)) {
        my $fastq_path = $fastq->name;
        if ($trim_paths) {
          $fastq_path =~ s/$root_trim// if $root_trim;
          $fastq_path = $ftp_root .'/'. $fastq_path if $ftp_root;
          $fastq_path =~ s{//+}{/}g;
        }
        push(@extra_header_lines, "\@CO\t$header_fastq_type = " . $fastq_path);
      }
    }
  }

}

my $reheader_object = ReseqTrack::Tools::ReheaderBam->new(
                  -input_files             => $input_file->name,
                  -working_dir             => $output_dir,
                  -header_lines_file       => $header_lines_file,
                  -samtools                => $samtools,
                  -job_name                => $name,
                  -extra_header_lines      => \@extra_header_lines,
                  );
$reheader_object->SQ_fields_hash('AS', $SQ_assembly);
$reheader_object->SQ_fields_hash('SP', $SQ_species);
$reheader_object->SQ_fields_hash('UR', $SQ_uri);

$reheader_object->options('reuse_old_header', $reuse_old_header ? 1 : 0);
$reheader_object->options('replace_PG', $replace_PG ? 1 : 0);
$reheader_object->options('replace_CO', $replace_CO ? 1 : 0);

$reheader_object->run;

$db->dbc->disconnect_when_inactive(0);
if($store){
  my $host = get_host_object($host_name, $db);

  my $bam = create_object_from_path($reheader_object->output_files->[0], $type_output, $host);
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

=pod

=head1 NAME

reseqtrack/scripts/process/run_reheader_bam.pl

=head1 SYNOPSIS

This script is uses samtools to reheader a bam file

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

  -name, name of a collection containing the bam file
  If name is a run_id or sample_id (or contains a run_id/sample_id), the run_meta_info table will be used to get some info

  -type_input, type of the collection containing the input bam file

  -type_output, collection type and file type when storing output bam file in the database

  -output_dir, base directory to hold files output bam files

  -header_lines_file, path to a file containing some / all of the header lines.
    if the -reuse_old_header is used, these lines will be appended to the header of the input bam
    if the -reuse_old_header is NOT used, these lines will become the new header

  -get_fastq_files, boolean.  When true, the collection table will be searched for fastq files with the correct run_id(s).
    These fastq files will be written to the @CO lines of the header.

  -type_fastq, type of fastq files when searching the collection table when the -get_fastq_files flag is used

  -samtools, path to the samtools executable.  The ReheaderBam class can guess at its location if it is not specified.

  -host_name, default is '1000genomes.ebi.ac.uk', needed for storing output bam file

  -store, boolean flag, to store output bam file in the database.

  -disable_md5, boolean flag, files written to the database will not have an md5

  -delete_input, boolean flag, to delete the original bam file (and index if it exists)
  and remove mark them as deleted in the database History table

  -directory_layout, specifies where the files will be located under output_dir.
      Tokens matching method names in RunMetaInfo will be substituted with that method's
      return value.

  -SQ_assembly, string e.g. 'GRCh37'. @SQ lines will be modified to contain the 'AS' field

  -SQ_species, string e.g. 'human'. @SQ lines will be modified to contain the 'SP' field

  -SQ_uri, string e.g. 'ftp://ftp.host/path/to/file'. @SQ lines will be modified to contain the 'UR' field

  -replace_PG, boolean.  When true, @PG lines in the input bam will not appear in the output bam

  -replace_CO, boolean.  When true, @CO lines in the input bam will not appear in the output bam

  -trim_paths, boolean.  Only applies when get_fastq_files is true.  File paths will be modified
      as described by -root_trim and -ftp_root

  -root_trim, string e.g. '/nfs/1000g-archive'.  When trim_paths is true, this string will be removed from fastq paths

  -ftp_root, string e.g. 'ftp://ftp.1000genomes.ebi.ac.uk'.  When trim_paths is true, this string will be added to fastq paths

  -run_id_regex, used to get run meta info and to assign fastq files as mate or frag.  Default is '[ESD]RR\d{6}'
  -study_id_regex, used to get run meta info.  Default is '[ESD]RS\d{6}'

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/process/run_bam_improvement.pl $DB_OPTS -name SRS000001
    -type_input BAM -type_output BAM_WITH_HEADER -output_dir /path/to/base_dir/
    -header_lines_file /path/to/header -get_fastq_files -type_fastq FILTERED_FASTQ
    -store -delete_input -directory_layout population/sample_id
    -SQ_assembly GRC37 -SQ_species human -replace_CO -replace_PG -reuse_old_header
    -trim_paths -root_trim /nfs/1000g-archive -ftp_root ftp://ftp.1000genomes.ebi.ac.uk

=cut

