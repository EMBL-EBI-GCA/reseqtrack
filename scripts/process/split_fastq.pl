#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(get_count_stats);
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);
use ReseqTrack::Tools::RunSplit;
use POSIX qw(ceil);
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $run_id;
my $type;
my $output_dir;
my $program_file;
my $max_reads;
my $help;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'run_id=s' => \$run_id,
  'type=s' => \$type,
  'max_reads=i' => \$max_reads,
  'output_dir=s' => \$output_dir,
  'program_file=s' => \$program_file,
  'help!'    => \$help,
    );

if ($help) {
    exec('perldoc', $0);
    exit(0);
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $ca = $db->get_CollectionAdaptor;

my $collection = $ca->fetch_by_name_and_type($run_id, $type);

throw("Failed to find a collection for ".$run_id." ".$type." from ".$dbname) 
    unless($collection);

my $rmi_a = $db->get_RunMetaInfoAdaptor;
my $run_meta_info = $rmi_a->fetch_by_run_id($run_id);

throw("Failed to find run_meta_info for $run_id from $dbname")
    unless($run_meta_info);

my $input_files = $collection->others;
my ($mate1, $mate2, $frag) = assign_files($input_files);

my %line_count_hash;

if ($mate1 && $mate2) {
  my ($read_count1,) = get_count_stats($mate1);
  my ($read_count2,) = get_count_stats($mate2);

  throw("read counts don't match: ".$mate1->name." $read_count1 "
        .$mate2->name." $read_count2")
        if ($read_count1 != $read_count2);

  my $num_output_files = ceil($read_count1 / $max_reads);
  my $line_count = ($num_output_files == 1) ? 0
                 : 4 * ceil($read_count1 / $num_output_files);
  $line_count_hash{$mate1->name} = $line_count;
  $line_count_hash{$mate2->name} = $line_count;
}
else {
  throw("Don't have mate pair for ".$mate1->name) if ($mate1);
  throw("Don't have mate pair for ".$mate2->name) if ($mate2);
}

if ($frag) {
  my ($read_count,) = get_count_stats($frag);
  my $num_output_files = ceil($read_count / $max_reads);
  my $line_count = ($num_output_files == 1) ? 0
                 : 4 * ceil($read_count / $num_output_files);
  $line_count_hash{$frag->name} = $line_count;
}
  

my $run_split = ReseqTrack::Tools::RunSplit->new(
    -program => $program_file,
    -working_dir => $output_dir,
    -line_count_hash => \%line_count_hash);

$run_split->run;

foreach my $file (@$input_files) {
  my $filename = $file->name;
  print "\ninput file:\n$filename:\n";
  my $outputs = $run_split->grouped_output_files($filename);
  print "output files:\n", join("\n", @$outputs), "\n";
}


=pod

=head1 NAME

ReseqTrack/scripts/process/split_fastq.pl

=head1 SYNOPSIS

    This script fetches a collection of fastq files. It splits the files into
    smaller files of similar sizes. The number of reads in the output files does
    not exceed a certain limit.

=head2 OPTIONS

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on
        -run_id, the run_id for the collection of fastq files
        -type, type of the collection, e.g. FILTERED_FASTQ
        -max_reads, integer, the number of reads in any output file will not exceed this value
        -output_dir, directory used for all output files
        -program_file, path to the split executable
        -help, flag to print this help and exit


=head1 Example:


$DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

    perl ReseqTrack/scripts/process/split_fastq.pl  $DB_OPTS -run_id ERR002097 -type FILTERED_FASTQ -max_reads 2000000
    -output_dir /path/to/dir -program_file ReseqTrack/c_code/split

=cut


