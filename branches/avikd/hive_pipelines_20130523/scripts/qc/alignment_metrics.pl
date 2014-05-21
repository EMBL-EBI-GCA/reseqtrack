#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::Tools::RunPicard;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::StatisticsUtils;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Loader::File;

use Getopt::Long;
use Data::Dumper;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;

my $java_path = '/usr/bin/java';
my $picard_dir;
my $name;
my $jvm_options = '-Xmx2g -Xms2g';
my $reference_sequence;
my $host_name         = '1000genomes.ebi.ac.uk';
my $keep_metrics_file = 0;
my $metrics_file_type = 'BAM_METRICS';
my $input_type;
my %options = (
  "validation_stringency" => "SILENT",
  "reference_sequence"    => $reference_sequence,
);

my @modes                   = ( 'alignment', 'rna' );
my $mode                    = $modes[0];
my $attribute_prefix_column = 'CATEGORY';

&GetOptions(
  'dbhost=s'                  => \$dbhost,
  'dbname=s'                  => \$dbname,
  'dbuser=s'                  => \$dbuser,
  'dbpass=s'                  => \$dbpass,
  'dbport=s'                  => \$dbport,
  'host=s'                    => \$host_name,
  'java_path=s'               => \$java_path,
  'picard_dir=s'              => \$picard_dir,
  'jvm_options=s'             => \$jvm_options,
  'reference_sequence=s'      => \$reference_sequence,
  'keep_metrics_file!'        => \$keep_metrics_file,
  'metrics_file_type=s'       => \$metrics_file_type,
  'input_type=s'              => \$input_type,
  'name=s'                    => \$name,
  'mode=s'                    => \$mode,
  'options=s'                 => \%options,
  'attribute_prefix_column=s' => \$attribute_prefix_column,
);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
);
$db->dbc->disconnect_when_inactive(1);
my $host = get_host_object( $host_name, $db );

my $collection_adaptor = $db->get_CollectionAdaptor;
my $collection =
  $collection_adaptor->fetch_by_name_and_type( $name, $input_type );

throw("Cannot find collection in db for $name $input_type") if ( !$collection );
throw("Collection $name $input_type must be of files")
  if ( $collection->table_name ne 'file' );
throw("Mode must be one of @modes , is $mode")
  if ( !grep { $mode eq $_ } @modes );

my $picard = ReseqTrack::Tools::RunPicard->new(
  -job_name     => $name,
  -program      => $java_path,
  -picard_dir   => $picard_dir,
  -jvm_options  => $jvm_options,
  -options      => \%options,
  -input_files  => [ $collection->others->[0]->name ],
  -keep_metrics => $keep_metrics_file,
);

my $picard_metrics;

if ( $mode eq 'alignment' ) {
  $picard_metrics = $picard->run_alignment_metrics;
}
elsif ( $mode eq 'rna' ) {
  $picard_metrics = $picard->run_rna_alignment_metrics;
}
else {
  throw("$mode did not match a run mode");
}

my $metrics_files = $picard->output_metrics_files;
my $statistics = [];

for my $metrics (@$picard_metrics) {
  my $prefix = $metrics->{$attribute_prefix_column};
  while ( my ( $key, $value ) = each %$metrics ) {
    next if $key eq $attribute_prefix_column;

    if ($prefix) {
      $key = join( '_', $prefix, $key );
    }

    push @$statistics, create_statistic_for_object( $collection, $key, $value )
      if ( defined $value );
  }
}

$collection->statistics($statistics);
$collection_adaptor->store_statistics($collection,1);

if ($keep_metrics_file) {
  create_output_records( $name, $metrics_file_type, $metrics_files );
}
else {
  for my $file (@$metrics_files) {
    delete_file($file);
  }
}

$picard->delete_files;

sub create_output_records {
  my ( $collection_name, $type, $file_names ) = @_;

  #create a File loader object for the files
  my $loader = ReseqTrack::Tools::Loader::File->new(
    -file            => $file_names,
    -do_md5          => 1,
    -hostname        => $host_name,
    -db              => $db,
    -assign_types    => 0,
    -check_types     => 0,
    -type            => $type,
    -update_existing => 1,
  );

  #Load the files into the database
  $loader->process_input();
  $loader->create_objects();
  $loader->sanity_check_objects();
  my $file_objects = $loader->load_objects();

  my $collection = ReseqTrack::Collection->new(
    -name       => $collection_name,
    -others     => $file_objects,
    -type       => $type,
    -table_name => 'file',
  );

  $db->get_CollectionAdaptor->store($collection);
}