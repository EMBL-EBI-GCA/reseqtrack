#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(get_count_stats create_object_from_path);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5);
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunSplit;
use POSIX qw(ceil);
use Getopt::Long;
use Data::Dumper;

$| = 1;

my $dbhost_r;
my $dbuser_r;
my $dbpass_r;
my $dbport_r = 4175;
my $dbname_r;
my $dbhost_w;
my $dbuser_w;
my $dbpass_w;
my $dbport_w = 4175;
my $dbname_w;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $run_id;
my $type_input;
my $type_output;
my $output_dir;
my $program_file;
my $max_reads;
my $help;

&GetOptions( 
  'dbhost_r=s'      => \$dbhost_r,
  'dbname_r=s'      => \$dbname_r,
  'dbuser_r=s'      => \$dbuser_r,
  'dbpass_r=s'      => \$dbpass_r,
  'dbport_r=s'      => \$dbport_r,
  'dbhost_w=s'      => \$dbhost_w,
  'dbname_w=s'      => \$dbname_w,
  'dbuser_w=s'      => \$dbuser_w,
  'dbpass_w=s'      => \$dbpass_w,
  'dbport_w=s'      => \$dbport_w,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'run_id=s' => \$run_id,
  'type_input=s' => \$type_input,
  'type_output=s' => \$type_output,
  'max_reads=i' => \$max_reads,
  'output_dir=s' => \$output_dir,
  'program_file=s' => \$program_file,
  'help!'    => \$help,
    );

if ($help) {
    exec('perldoc', $0);
    exit(0);
}

my $db_r = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost_r,
  -user   => $dbuser_r,
  -port   => $dbport_r,
  -dbname => $dbname_r,
  -pass   => $dbpass_r,
    );

my $db_w = $store ? ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost_w,
  -user   => $dbuser_w,
  -port   => $dbport_w,
  -dbname => $dbname_w,
  -pass   => $dbpass_w,
    ) : undef;

my $ca_r = $db_r->get_CollectionAdaptor;

my $collection = $ca_r->fetch_by_name_and_type($run_id, $type_input);

throw("Failed to find a collection for ".$run_id." ".$type_input." from ".$dbname_r) 
    if(!$collection);

my $rmi_a = $db_r->get_RunMetaInfoAdaptor;
my $run_meta_info = $rmi_a->fetch_by_run_id($run_id);

throw("Failed to find run_meta_info for $run_id from $dbname_r")
    if (!$run_meta_info);

my $input_files = $collection->others;
my ($mate1, $mate2, $frag) = assign_files($input_files);
my ($mate1_path, $mate2_path, $frag_path) = map {$_ ? $_->name : ''} ($mate1, $mate2, $frag);

my %line_count_hash;

if ($mate1 && $mate2) {
  my ($read_count1,) = get_count_stats($mate1);
  my ($read_count2,) = get_count_stats($mate2);

  throw("read counts don't match: ".$mate1_path." $read_count1 "
        .$mate2_path." $read_count2")
        if ($read_count1 != $read_count2);

  my $num_output_files = ceil($read_count1 / $max_reads);
  my $line_count = ($num_output_files == 1) ? 0
                 : 4 * ceil($read_count1 / $num_output_files);
  $line_count_hash{$mate1_path} = $line_count;
  $line_count_hash{$mate2_path} = $line_count;
}
else {
  throw("Don't have mate pair for ".$mate1_path) if ($mate1);
  throw("Don't have mate pair for ".$mate2_path) if ($mate2);
}

if ($frag) {
  my ($read_count,) = get_count_stats($frag);
  my $num_output_files = ceil($read_count / $max_reads);
  my $line_count = ($num_output_files == 1) ? 0
                 : 4 * ceil($read_count / $num_output_files);
  $line_count_hash{$frag_path} = $line_count;
}
  

my $run_split = ReseqTrack::Tools::RunSplit->new(
    -program => $program_file,
    -working_dir => $output_dir,
    -line_count_hash => \%line_count_hash);

$run_split->run;

if ($store) {
  my $host = get_host_object($host_name, $db_w);
  my $ca_w = $db_w->get_CollectionAdaptor;

  my $output_file_hash = $run_split->output_file_hash;
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
        my $md5 = run_md5($file_path);
        $file->md5($md5);
        push(@files, $file);
      }
      my $collection_name = $run_id . '_m' . $label;
      my $collection = ReseqTrack::Collection->new(
                    -name => $collection_name, -type => $type_output,
                    -table_name => 'file', -others => \@files);
      push(@split_collections, $collection);
    }
  }

  if ($frag) {
    while (my ($label, $file_path) = each(%{$output_file_hash->{$frag_path}})) {
      my $file= create_object_from_path($file_path, $type_output, $host);
      my $md5 = run_md5($file_path);
      $file->md5($md5);
      my $collection_name = $run_id . '_f' . $label;
      my $collection = ReseqTrack::Collection->new(
                    -name => $collection_name, -type => $type_output,
                    -table_name => 'file', -others => $file);
      push(@split_collections, $collection);
    }
  }

      my $collection = ReseqTrack::Collection->new(
                    -name => $run_id, -type => $type_output,
                    -table_name => 'collection', -others => \@split_collections);
      $ca_w->store($collection);
}



=pod

=head1 NAME

ReseqTrack/scripts/process/split_fastq.pl

=head1 SYNOPSIS

    This script fetches a collection of fastq files. It splits the files into
    smaller files of similar sizes. The number of reads in the output files does
    not exceed a certain limit.

=head2 OPTIONS

      options for the database that stores in input fastq files:

        -dbhost_r, the name of the mysql-host
        -dbname_r, the name of the mysql database
        -dbuser_r, the name of the mysql user
        -dbpass_r, the database password if appropriate
        -dbport_r, the port the mysql instance is running on

      options for the database to hold the output fastq files:

        -dbhost_w, the name of the mysql-host
        -dbname_w, the name of the mysql database
        -dbuser_w, the name of the mysql user
        -dbpass_w, the database password if appropriate
        -dbport_w, the port the mysql instance is running on
        -store, flag to store output fastq files to the database

      other options:

        -run_id, the run_id for the collection of fastq files
        -type_input, type of the collection of input files, e.g. FILTERED_FASTQ
        -type_output, type of the collection of output files, e.g. SPLIT_FASTQ
        -max_reads, integer, the number of reads in any output file will not exceed this value
        -output_dir, directory used for all output files
        -program_file, path to the split executable
        -host_name, name of the host object to associate with the output files
        -help, flag to print this help and exit


=head1 Example:


    $DB_R_OPTS="-dbhost_r mysql-host -dbuser_r rw_user -dbpass_r **** -dbport_r 4197 -dbname_r my_database_1"
    $DB_W_OPTS="-dbhost_w mysql-host -dbuser_w rw_user -dbpass_w **** -dbport_w 4197 -dbname_w my_database_2"

    perl ReseqTrack/scripts/process/split_fastq.pl  $DB_R_OPTS $DB_W_OPTS
      -run_id ERR002097 -type_input FILTERED_FASTQ -type_output SPLIT_FASTQ -max_reads 2000000
      -output_dir /path/to/dir -program_file ReseqTrack/c_code/split -store

=cut


