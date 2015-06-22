#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(get_count_stats create_object_from_path);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5);
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw( create_directory_path );
use ReseqTrack::Tools::RunSplit;
use POSIX qw(ceil);
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $host_name = '1000genomes.ebi.ac.uk';
my $type_input;
my $type_output;
my $type_collection;
my $max_reads;
my $max_bases;
my $help;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'host_name=s' => \$host_name,
  'type_input=s' => \$type_input,
  'type_output=s' => \$type_output,
  'type_collection=s' => \$type_collection,
  'max_reads=i' => \$max_reads,
  'max_bases=i' => \$max_bases,
  'help!'    => \$help,
    );

if ($help) {
    exec('perldoc', $0);
    exit(0);
}

throw("specify only one of max_reads or max_bases")
  if ($max_reads && $max_bases);

throw("must specify type of output")
  if (!$type_output);
throw("must specify type of collection")
  if (!$type_collection);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );
$db->dbc->disconnect_when_inactive(1);


my $ca = $db->get_CollectionAdaptor;

my %existing_names;
foreach my $collection (@{$ca->fetch_by_type($type_collection)}) {
    $existing_names{$collection->name} = 1;
}

my @input_collections = grep {! $existing_names{$_->name}} @{$ca->fetch_by_type($type_input)};

foreach my $input_collection (@input_collections) {

  my $input_files = $input_collection->others;
  my ($mate1, $mate2, $frag) = assign_files($input_files);
  my ($mate1_path, $mate2_path, $frag_path) = map {$_ ? $_->name : ''} ($mate1, $mate2, $frag);


  my @split_collections;
  if ($mate1 || $mate2) {
    throw("Don't have mate pair for ".$mate1_path) if (!$mate2);
    throw("Don't have mate pair for ".$mate2_path) if (!$mate1);

    my ($read_count1, $base_count1) = get_count_stats($mate1);
    my ($read_count2, $base_count2) = get_count_stats($mate2);

    throw("read counts don't match: ".$mate1_path." $read_count1 "
          .$mate2_path." $read_count2")
          if ($read_count1 != $read_count2);

    my $num_output_files = $max_reads ? ceil($read_count1 / $max_reads)
                        : ceil( 0.5 * ($base_count1 +$base_count2) / $max_bases );

    my @collection_names;
    if ($num_output_files == 1) {
        push(@collection_names, $input_collection->name . '_MATE');
    }
    else {
      my $reads_per_split = ($num_output_files == 1) ? 0
                       : ceil($read_count1 / $num_output_files);

      foreach my $i (0..$num_output_files-1) {
        my $read_first = $i * $reads_per_split + 1;
        my $read_last = $i == $num_output_files -1 
                      ? $read_count1
                      : $read_first + $reads_per_split -1;
        push(@collection_names, $input_collection->name . "_MATE$read_first-$read_last");
      }
    }
    foreach my $collection_name (@collection_names) {
      my $collection = ReseqTrack::Collection->new(
                    -name => $collection_name, -type => $type_output,
                    -table_name => 'file', -others => [$mate1, $mate2]);
      push(@split_collections, $collection);
    }

  }

  if ($frag) {
    my ($read_count, $base_count) = get_count_stats($frag);

    my $num_output_files = $max_reads ? ceil($read_count / $max_reads)
                        : ceil( $base_count / $max_bases );

    my @collection_names;
    if ($num_output_files == 1) {
        push(@collection_names, $input_collection->name . '_FRAG');
    }
    else {
      my $reads_per_split = ($num_output_files == 1) ? 0
                       : ceil($read_count / $num_output_files);

      foreach my $i (0..$num_output_files-1) {
        my $read_first = $i * $reads_per_split + 1;
        my $read_last = $i == $num_output_files -1 
                      ? $read_count
                      : $read_first + $reads_per_split -1;
        push(@collection_names, $input_collection->name . "_FRAG$read_first-$read_last");
      }
    }
    foreach my $collection_name (@collection_names) {
      my $collection = ReseqTrack::Collection->new(
                    -name => $collection_name, -type => $type_output,
                    -table_name => 'file', -others => [$frag]);
      push(@split_collections, $collection);
    }

  }

  my $collection_of_collections = ReseqTrack::Collection->new(
                -name => $input_collection->name, -type => $type_collection,
                -table_name => 'collection', -others => \@split_collections);
  $ca->store($collection_of_collections);
}




=pod

=head1 NAME

ReseqTrack/scripts/process/make_split_fastq_collections.pl

=head1 SYNOPSIS

    This script fetches collections of fastq files. It calculates the best way to split
    the files for alignment. It writes new collections to the databases with names like
    NAME_MATE1-1000 and NAME_FRAG1-1000

=head2 OPTIONS

      database options:

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on
        -host_name, name of the host object to associate with the output files
            (default is 1000genomes.ebi.ac.uk)

      other options:

        -type_input, type of the input collection of files, e.g. FILTERED_FASTQ
        -type_output, type of the output collection of files
            (many of these are created for each input collection)
        -type_collection, type of the collection of collections of files, e.g. FASTQ_CHUNK_SET
            (one of these is created for each input collection)

        -max_reads, integer, the number of reads will not exceed this value
        -max_bases, integer, the number of bases will not exceed this value
            Specify only one of -max_reads or -max_bases

        -help, flag to print this help and exit


=head1 Example:


    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

    perl ReseqTrack/scripts/process/split_fastq.pl  $DB_OPTS
      -type_input FILTERED_FASTQ -type_output FASTQ_CHUNK -type_collection FASTQ_CHUNK_SET -max_reads 2000000

=cut


