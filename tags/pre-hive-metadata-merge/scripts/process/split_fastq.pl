#!/usr/bin/env perl

use strict;
use warnings;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(get_count_stats create_object_from_path);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5 check_directory_exists);
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw( create_directory_path );
use ReseqTrack::Tools::RunSplit;
use File::Basename qw(fileparse);
use File::Copy qw (move copy);
use POSIX qw(ceil);
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $disable_md5;
my $run_id;
my $type_input;
my $type_output;
my $type_collection;
my $output_dir;
my $program_file;
my $max_reads;
my $max_bases;
my $help;
my $directory_layout = 'population/sample_id/run_id';
my $no_split_strategy = 'copy';
my $run_id_regex = '[ESD]RR\d{6}';

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'disable_md5!' => \$disable_md5,
  'run_id=s' => \$run_id,
  'type_input=s' => \$type_input,
  'type_output=s' => \$type_output,
  'type_collection=s' => \$type_collection,
  'max_reads=i' => \$max_reads,
  'max_bases=i' => \$max_bases,
  'output_dir=s' => \$output_dir,
  'program_file=s' => \$program_file,
  'help!'    => \$help,
  'directory_layout=s' => \$directory_layout,
  'no_split_strategy=s' => \$no_split_strategy,
  'run_id_regex=s' => \$run_id_regex,
    );

if ($help) {
    exec('perldoc', $0);
    exit(0);
}

throw("specify only one of max_reads or max_bases")
  if ($max_reads && $max_bases);

throw("no_split_strategy must be move, copy, link or ignore")
  if (!grep {$no_split_strategy eq $_} qw(move copy link ignore));

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );
$db->dbc->disconnect_when_inactive(1);

my $ca = $db->get_CollectionAdaptor;

my $collection = $ca->fetch_by_name_and_type($run_id, $type_input);

throw("Failed to find a collection for ".$run_id." ".$type_input." from ".$dbname) 
    if(!$collection);

my $rmi_a = $db->get_RunMetaInfoAdaptor;
my $run_meta_info = $rmi_a->fetch_by_run_id($run_id);

throw("Failed to find run_meta_info for $run_id from $dbname")
    if (!$run_meta_info);

my $input_files = $collection->others;
my @regexs = (qr/$run_id_regex\S*_1\.(\w+\.)*fastq(\.gz)?$/i,
              qr/$run_id_regex\S*_2\.(\w+\.)*fastq(\.gz)?$/i,
              qr/$run_id_regex\S*\.fastq(\.gz)?$/i);
my ($mate1, $mate2, $frag) = assign_files($input_files, \@regexs);
my ($mate1_path, $mate2_path, $frag_path) = map {$_ ? $_->name : ''} ($mate1, $mate2, $frag);

my %line_count_hash;
my @files_no_split;

if ($mate1 && $mate2) {
  my ($read_count1, $base_count1) = get_count_stats($mate1);
  my ($read_count2, $base_count2) = get_count_stats($mate2);

  throw("read counts don't match: ".$mate1_path." $read_count1 "
        .$mate2_path." $read_count2")
        if ($read_count1 != $read_count2);

  my $num_output_files = $max_reads ? ceil($read_count1 / $max_reads)
                      : ceil( 0.5 * ($base_count1 +$base_count2) / $max_bases );
  if ($num_output_files ==1) {
    push(@files_no_split, $mate1_path, $mate2_path);
  }
  else {
    my $line_count = 4 * ceil($read_count1 / $num_output_files);
    $line_count_hash{$mate1_path} = $line_count;
    $line_count_hash{$mate2_path} = $line_count;
  }
}
else {
  throw("Don't have mate pair for ".$mate1_path) if ($mate1);
  throw("Don't have mate pair for ".$mate2_path) if ($mate2);
}

if ($frag) {
  my ($read_count, $base_count) = get_count_stats($frag);
  my $num_output_files = $max_reads ? ceil($read_count / $max_reads)
                      : ceil($base_count / $max_bases);
  if ($num_output_files ==1) {
    push(@files_no_split, $frag_path);
  }
  else {
    my $line_count = 4 * ceil($read_count / $num_output_files);
    $line_count_hash{$frag_path} = $line_count;
  }
}

my $full_output_dir = create_directory_path($run_meta_info, $directory_layout, $output_dir);

my $output_file_hash;
if (scalar keys %line_count_hash) {
  my $run_split = ReseqTrack::Tools::RunSplit->new(
      -program => $program_file,
      -working_dir => $full_output_dir,
      -line_count_hash => \%line_count_hash,
      -job_name => $run_id);

  $run_split->run;
  $output_file_hash = $run_split->output_file_hash;
}

if ($no_split_strategy ne 'ignore' && scalar @files_no_split) {
  check_directory_exists($full_output_dir);
  foreach my $filepath (@files_no_split) {
    my ($basename, $input_dir, $suffix) = fileparse($filepath, qr/\.[^\/]+/);
    my $output_file = "$full_output_dir/$basename.not_split$suffix";
    $output_file =~ s{//}{/}g;
    $output_file_hash->{$filepath}->{''} = $output_file;
    if ($no_split_strategy eq 'move') {
      move ($filepath, $output_file)
          or throw "move failed: $!";
    }
    elsif ($no_split_strategy eq 'copy') {
      copy ($filepath, $output_file)
          or throw "copy failed: $!";
    }
    elsif ($no_split_strategy eq 'link') {
      symlink ($filepath, $output_file)
          or throw "link failed: $!";
    }
  }
}

if ($store) {
  my $host = get_host_object($host_name, $db);

  my @split_collections;

  if ($mate1 && $mate2) {
    my @labels_mate1 = keys %{$output_file_hash->{$mate1_path}};
    throw("different number of output files for mate1 and mate2")
      if (@labels_mate1 != keys %{$output_file_hash->{$mate2_path}});
    foreach my $label (@labels_mate1) {
      throw("no matching mate2 file with label $label") if (!$output_file_hash->{$mate2_path}->{$label});
      my @files;
      foreach my $mate ($mate1_path, $mate2_path){
        my $file_path = $output_file_hash->{$mate}->{$label};
        my $file = create_object_from_path($file_path, $type_output, $host);
        if (! $disable_md5) {
          $file->md5(run_md5($file_path));
        }
        my $statistic = ReseqTrack::Statistic->new(
                    -attribute_name => 'paired_length',
                    -attribute_value => $run_meta_info->paired_length);
        $file->statistics($statistic);
        push(@files, $file);
      }
      my $collection_name = $run_id . '_m';
      $collection_name .= $label if $label;
      my $collection = ReseqTrack::Collection->new(
                    -name => $collection_name, -type => $type_output,
                    -table_name => 'file', -others => \@files);

      push(@split_collections, $collection);
    }
  }

  if ($frag) {
    while (my ($label, $file_path) = each(%{$output_file_hash->{$frag_path}})) {
      my $file= create_object_from_path($file_path, $type_output, $host);
      if (! $disable_md5) {
        $file->md5(run_md5($file_path));
      }
      my $collection_name = $run_id . '_f';
      $collection_name .= $label if $label;
      my $collection = ReseqTrack::Collection->new(
                    -name => $collection_name, -type => $type_output,
                    -table_name => 'file', -others => $file);
      push(@split_collections, $collection);
    }
  }

  my $collection = ReseqTrack::Collection->new(
                -name => $run_id, -type => $type_collection,
                -table_name => 'collection', -others => \@split_collections);
  $ca->store($collection);
}



=pod

=head1 NAME

ReseqTrack/scripts/process/split_fastq.pl

=head1 SYNOPSIS

    This script fetches a collection of fastq files. It splits the files into
    smaller files of similar sizes. The number of reads or number of bases in the
    output files does not exceed a certain limit.

=head2 OPTIONS

      database options:

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on
        -host_name, name of the host object to associate with the output files
            (default is 1000genomes.ebi.ac.uk)
        -store, flag to store output fastq files to the database
        -disable_md5, boolean flag, files written to the database will not have an md5

      other options:

        -run_id, the run_id for the collection of fastq files
        -type_input, type of the collection of input files, e.g. FILTERED_FASTQ
        -type_output, type of the output files and collection of output files, e.g. FASTQ_CHUNK
        -type_collection, type of the collection of collections of output files, e.g. FASTQ_CHUNK_SET

        -max_reads, integer, the number of reads in any output file will not exceed this value
        -max_bases, integer, the number of bases in any output file will not exceed this value
            Specify only one of -max_reads or -max_bases

        -output_dir, base directory used for all output files
        -directory_layout, specifies where the files will be located under output_dir.
              Tokens matching method names in RunMetaInfo will be substituted with that method's return value.
              Default value is population/sample_id/run_id
        -program_file, path to the split executable

        -no_split_strategy, must be either 'move', 'copy', 'link' or 'ignore'
              Tells the script what to do with the original fastq file if it does not need to be split

        -help, flag to print this help and exit


=head1 Example:


    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

    perl ReseqTrack/scripts/process/split_fastq.pl  $DB_OPTS
      -run_id ERR002097 -type_input FILTERED_FASTQ -type_output FASTQ_CHUNK -type_collection FASTQ_CHUNK_SET
      -max_reads 2000000 -output_dir /path/to/dir -program_file ReseqTrack/c_code/split -store
      -no_split_strategy link

=cut


