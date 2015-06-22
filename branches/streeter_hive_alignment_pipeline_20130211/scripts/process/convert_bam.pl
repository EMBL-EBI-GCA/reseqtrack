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
use ReseqTrack::Tools::GeneralUtils qw(execute_system_command);
use Getopt::Long;

my $dbhost;
my $dbname;
my $dbuser;
my $dbpass;
my $dbport;
my $name;
my $type_input;
my $type_output;
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
my $gzip_output;
my $help;

&GetOptions(
  'dbhost=s'           => \$dbhost,
  'dbname=s'           => \$dbname,
  'dbuser=s'           => \$dbuser,
  'dbpass=s'           => \$dbpass,
  'dbport=s'           => \$dbport,
  'name=s'             => \$name,
  'type_input=s'       => \$type_input,
  'type_output=s'      => \$type_output,
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
  'gzip_output!'       => \$gzip_output,
  'help!'              => \$help,
);

if ($help) {
  exec( 'perldoc', $0 );
  exit(0);
}

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

  if ($gzip_output) {
    for ( my $i = 0 ; $i < scalar(@file_paths) ; $i++ ) {
      execute_system_command( 'gzip ' . $file_paths[$i] );
      $file_paths[$i] .= '.gz';
    }
  }

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

=pod

=head1 NAME

ReseqTrack/scripts/process/convert_bam.pl

=head1 SYNOPSIS

    This script fetches a collection containing a single BAM file. It converts the BAM file into a different format
    and creates a collection for that file.

=head2 OPTIONS

      database options:

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on
        -host_name, name of the host object to associate with the output files
            (default is 1000genomes.ebi.ac.uk)
        -store, flag to store output files to the database
        -disable_md5, boolean flag, files written to the database will not have an md5

      general options:
        -name, name of the collection to load, also used as the name of the output collection
        -type_input, type of the input collection
        -type_output, type for the output file and collection
        -output_dir, root location for the output file
        -directory_layout, how to place the output below the root directory. Can use run or sample information.
        -chromosome_file, location of a 2 column, tab-separated file listing the chromosomes and their lengths. Required for producing bedgraph, bigwig and bigbed output;
        -run_id_regex, how to identify a name as a run_id, defaults to '[ESD]RR\d{6}'
        -sample_id_regex, how to identify a name as a sample_id, defaults to '[ESD]RS\d{6}'
        -output_format, which conversion to run bg for bedgraph, bw for BigWig, bed for bed and bb for bigbed;
        -options hash of options to pass to the conversion;
        -gzip_output flag to gzip the output. This is not sensible for the compressed formats BigBed and BigWig, but the script allows it
        -help, flag to print this help and exit
        
      options hash:
        -options split_reads=1, will cause bedTools executables to be run with -split, as appropriate for RNA-seq reads

=head1 Dependencies:
  
    The bedTools suite executables should be in the PATH, genomeCoverageBed and bamToBed are used.
    Jim Kent's bedGraphToBigWig and bedToBigBed should also be in the PATH.

=head1 Example:


    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

    perl ReseqTrack/scripts/process/convert_bam.pl  $DB_OPTS
      -type_input BAM -type_output BIG_WIG -output_format bw
      -output_dir /path/to/dir -chromosome_file /path/to/file
      -store


=cut
