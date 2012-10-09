use strict;
use warnings;

use ReseqTrack::Tools::ConvertBam;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use ReseqTrack::Tools::Loader::File;
use Getopt::Long;

my $dbhost;
my $dbname;
my $dbuser;
my $dbpass;
my $dbport;
my $name;
my $type_input;
my $type_output;
my $type_index;
my $output_dir;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $disable_md5;
my $directory_layout;
my $chromosome_file;
my $run_id_regex    = '[ESD]RR\d{6}';
my $sample_id_regex = '[ESD]RS\d{6}';
my $output_format;
my %options;

&GetOptions(
  'dbhost=s'           => \$dbhost,
  'dbname=s'           => \$dbname,
  'dbuser=s'           => \$dbuser,
  'dbpass=s'           => \$dbpass,
  'dbport=s'           => \$dbport,
  'name=s'             => \$name,
  'type_input=s'       => \$type_input,
  'type_output=s'      => \$type_output,
  'type_index=s'       => \$type_index,
  'output_dir=s'       => \$output_dir,
  'host_name=s'        => \$host_name,
  'store!'             => \$store,
  'disable_md5!'       => \$disable_md5,
  'directory_layout=s' => \$directory_layout,
  'chromosome_file=s'  => \$chromosome_file,
  'options=s'          => \%options,
  'run_id_regex=s'     => \$run_id_regex,
  'sample_id_regex=s'  => \$sample_id_regex,
  'output_format=s'    => \$output_format,
);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
);
$db->dbc->disconnect_when_inactive(1);
my $ca = $db->get_CollectionAdaptor;
my $fa = $db->get_FileAdaptor;

my $collection = $ca->fetch_by_name_and_type( $name, $type_input );
throw("Failed to find a collection for "
    . $name . " "
    . $type_input
    . " from "
    . $dbname )
  unless ($collection);

my $input_files = $collection->others;
my @input_filepaths = map { $_->{'name'} } @$input_files;

if ($directory_layout) {
  my $rmia = $db->get_RunMetaInfoAdaptor;
  my $run_meta_info;

  if ( $name =~ /$run_id_regex/ ) {
    $run_meta_info = $rmia->fetch_by_run_id($&);
  }
  elsif ( $name =~ /$sample_id_regex/ ) {
    my $rmi_list = $rmia->fetch_by_sample_id($&);
    $run_meta_info = $rmi_list->[0] if (@$rmi_list);
  }

  if ($run_meta_info) {
    $output_dir =
      create_directory_path( $run_meta_info, $directory_layout, $output_dir );
  }
}

my $converter = ReseqTrack::Tools::ConvertBam->new(
  -input_files     => \@input_filepaths,
  -working_dir     => $output_dir,
  -job_name        => $name,
  -options         => \%options,
  -chromosome_file => $chromosome_file,
  -output_format   => $output_format,
  -options         => \%options,
);

$converter->run;

if ($store) {
  my $host = get_host_object( $host_name, $db );

  my @file_paths = @{ $converter->output_files };

  if (@file_paths) {
    foreach my $path (@file_paths) {
      throw("database already has file with name $path")
        if ( $fa->fetch_by_name($path) );
    }
    my $files =
      create_objects_from_path_list( \@file_paths, $type_output, $host );
    if ( !$disable_md5 ) {
      foreach my $file (@$files) {
        $file->md5( run_md5( $file->name ) );
      }
    }

    my $collection = ReseqTrack::Collection->new(
      -name   => $name,
      -type   => $type_output,
      -others => $files
    );

    $ca->store($collection);
  }

}

