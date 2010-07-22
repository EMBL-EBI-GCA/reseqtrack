#!/usr/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::RunAlignment::BFAST;
use Getopt::Long;
use Data::Dumper;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $run_id;
my $type;
my $output_dir;
my $build;
my $read_length;
my $working_dir;
&GetOptions(
             'dbhost=s'      => \$dbhost,
             'dbname=s'      => \$dbname,
             'dbuser=s'      => \$dbuser,
             'dbpass=s'      => \$dbpass,
             'dbport=s'      => \$dbport,
             'run_id=s'      => \$run_id,
             'type=s'        => \$type,
             'output_dir=s'  => \$output_dir,
             'build=s'       => \$build,
             'read_length=s' => \$read_length,
);

my $db =
  ReseqTrack::DBSQL::DBAdaptor->new(
                                     -host   => $dbhost,
                                     -user   => $dbuser,
                                     -port   => $dbport,
                                     -dbname => $dbname,
                                     -pass   => $dbpass,
  );

my $ca = $db->get_CollectionAdaptor;
my $collection = $ca->fetch_by_name_and_type( $run_id, $type );

throw(   "Failed to find a collection for  $run_id   $type  from     . $dbname \n")
  unless ($collection);

my $bfast = "/nfs/1000g-work/G1K/work/bin/bfast-0.6.4d/bin/bfast";
my $preprocess_exe =
"/nfs/1000g-work/G1K/work/bin/bfast-0.6.4d/bin/1kg2bfastfastq/1kg_2_bfast_fastq";
my $reference     = "";
my $short_indexes = "/nfs/1000g-work/G1K/work/reference/BFAST/short/";
my $long_indexes  = "/nfs/1000g-work/G1K/work/reference/BFAST/long/";
my $short_ref     = "${short_indexes}human_g1k_v37.fa";
my $long_ref      = "${long_indexes}human_g1k_v37.fa";
my $samtools      = "/nfs/1000g-work/G1K/work/bin/samtools/samtools";

my $run_alignment = ReseqTrack::Tools::RunAlignment::BFAST->new(
               -input           => $collection,
               -name            => $run_id,
               -read_length     => $read_length,
               -working_dir     => '/home/smithre/RUN_FAST/grc',
               -program         => $bfast,
               -preprocess_exe  => $preprocess_exe,
               -reference       => $reference,
               -short_reference => $short_ref,
               -long_reference  => $long_ref,
               -samtools        => $samtools,
);

$run_alignment->run();


