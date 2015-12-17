#!/usr/bin/env perl

use strict;
use warnings;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list);
use ReseqTrack::Tools::RunWiggler;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::MetaDataUtils
  qw(create_directory_path fetch_metadata_object);
use ReseqTrack::Tools::GeneralUtils
  qw(execute_system_command execute_pipe_system_command);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5);
use ReseqTrack::Tools::AttributeUtils qw(create_attribute_for_object);
use ReseqTrack::Tools::Loader::File;

use Statistics::Descriptive;

use File::Copy "move";
use File::Basename qw(fileparse dirname);

use Data::Dumper;
use Getopt::Long;

$| = 1;

my ( $dbhost, $dbuser, $dbpass, $dbport, $dbname );
my ( $name, $type_input, $type_output, $output_dir, $directory_layout,
    $gzip_output );
my $program;
my $verbose                 = 0;
my $host_name               = '1000genomes.ebi.ac.uk';
my $run_id_regex            = '[ESD]RR\d{6}';
my $sample_id_regex         = '[ESD]RS\d{6}';
my $preserve_paths          = 0;
my $fragment_size_stat_name = '';
my $dedupe                  = 0;
my $store                   = 0;
my $disable_md5             = 0;
my $samtools_path           = 0;
my $chrom_sizes_file;
my $mcr_root;

my %options;

my %config;
my $bedGraph_to_bigWig_path;

&GetOptions(
    'dbhost=s'                  => \$dbhost,
    'dbname=s'                  => \$dbname,
    'dbuser=s'                  => \$dbuser,
    'dbpass=s'                  => \$dbpass,
    'dbport=s'                  => \$dbport,
    'name=s'                    => \$name,
    'type_input=s'              => \$type_input,
    'type_output=s'             => \$type_output,
    'output_dir=s'              => \$output_dir,
    'directory_layout=s'        => \$directory_layout,
    'gzip_output!'              => \$gzip_output,
    'store!'                    => \$store,
    'program=s'                 => \$program,
    'disable_md5!'              => \$disable_md5,
    'options=s'                 => \%options,
    'preserve_paths'            => \$preserve_paths,
    'fragment_size_stat_name=s' => \$fragment_size_stat_name,
    'bedGraphToBigWig_path=s'   => \$bedGraph_to_bigWig_path,
    'chrom_sizes_file=s'        => \$chrom_sizes_file,
    'dedupe!'                   => \$dedupe,
    'host=s'                    => \$host_name,
    'mcr_root=s'                => \$mcr_root,
    'samtools_path=s'           => \$samtools_path,
) or die($!);

throw("Must specify an output directory ")
  if ( !$output_dir && !$preserve_paths );
throw("Must specify an output type") if ( !$type_output );

# dba + adaptors
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);
throw("Could not create db adaptor") if ( !$db );

$db->dbc->disconnect_when_inactive(1);
my $ca   = $db->get_CollectionAdaptor;
my $fa   = $db->get_FileAdaptor;
my $rmia = $db->get_RunMetaInfoAdaptor;

#inputs
my $collection = $ca->fetch_by_name_and_type( $name, $type_input );
throw(  "Failed to find a collection for "
      . $name . " "
      . $type_input
      . " from "
      . $dbname )
  unless ($collection);

my $input_files     = $collection->others;
my $job_name        = $name;
my @input_filepaths = map { $_->{'name'} } @$input_files;

#
my $mcr_jre      = "$mcr_root/sys/java/jre/glnxa64/jre/lib/amd64";
my $xapplres_dir = "$mcr_root/X11/app-defaults";
my $ld_lib_path  = join( ':',
    $ENV{LD_LIBRARY_PATH},     "$mcr_root/runtime/glnxa64",
    "$mcr_root/bin/glnxa64",   "$mcr_root/sys/os/glnxa64",
    "$mcr_jre/native_threads", "$mcr_jre/server",
    $mcr_jre );

$ENV{MCRROOT}         = $mcr_root;
$ENV{MCRJRE}          = $mcr_jre;
$ENV{XAPPLRESDIR}     = $xapplres_dir;
$ENV{LD_LIBRARY_PATH} = $ld_lib_path;

if ($samtools_path){
  $ENV{PATH} = dirname($samtools_path).':'.$ENV{PATH};
}

my $meta_info = fetch_metadata_object( $name, $db );
if ($preserve_paths) {
    my ( $filename, $basepath, $suffix ) = fileparse( $input_filepaths[0] );
    print Dumper( [ $filename, $basepath, $suffix ] );
    $output_dir = $basepath;
    $job_name   = $filename;
    $job_name =~ s/\.bam$//;
}
else {
    $output_dir =
      create_directory_path( $meta_info, $directory_layout, $output_dir )
      if ($directory_layout);
}

print 'Found input files ' . join( ', ', @input_filepaths ) . $/ if $verbose;

if ($fragment_size_stat_name) {
    my $fragment_length;

    # from stats on BAM file

    my $stats = $collection->attributes();

    for my $stat (@$stats) {
        if ( $stat->attribute_name eq $fragment_size_stat_name ) {

# ppqt can find multiple peaks. the first is the most likely, so we use that one
            ($fragment_length) = grep { $_ >= 1 } split /,/,
              $stat->attribute_value();
            if ( !defined $fragment_length ) {
                $fragment_length = $stat->attribute_value;
            }

        }
    }

    throw(
"Could not find statistic $fragment_size_stat_name for collection $name $type_input"
    ) unless defined $fragment_length;

    if ( $fragment_length < 1 ) {
        $fragment_length = 1;    #minimum possible value
    }
    print STDERR "Effective fragment length $fragment_length $/";

    $options{l} = $fragment_length;
}

$config{-program}       = $program;
$config{-input_files}   = \@input_filepaths;
$config{-echo_cmd_line} = $verbose;
$config{-working_dir}   = $output_dir;
$config{-job_name}      = $job_name;
$config{-options}       = \%options;
$config{-dedupe}        = $dedupe;
$config{-samtools_path} = $samtools_path;

if ($bedGraph_to_bigWig_path) {
    $config{-output_format}         = 'bw';
    $config{-bedgraphtobigwig_path} = $bedGraph_to_bigWig_path;
    $config{-chrom_sizes_file}      = $chrom_sizes_file;
}
else {
    $config{-output_format} = 'bg';
}

# run analysis
my $wiggler = ReseqTrack::Tools::RunWiggler->new(%config);

$wiggler->run();
my @files = @{ $wiggler->output_files };

if ($gzip_output) {
    for ( my $i = 0 ; $i < scalar(@files) ; $i++ ) {
        execute_system_command("gzip -f $files[$i]");
        $files[$i] .= '.gz';
    }
}

for my $output (@files) {
    print "OUTPUT:" . $output . $/;
    throw("Cannot find output file: $output") if ( !-e $output );
}

create_output_records( $name, $type_output, \@files );

sub get_run_meta_info {
    my ( $name, $rmia, $run_id_regex, $sample_id_regex ) = @_;

    my $run_meta_info;
    if ( $name =~ /$run_id_regex/ ) {
        $run_meta_info = $rmia->fetch_by_run_id($name);
    }
    elsif ( $name =~ /$sample_id_regex/ ) {
        my $rmi_list = $rmia->fetch_by_sample_id($name);
        $run_meta_info = $rmi_list->[0] if (@$rmi_list);
    }
    else {
        throw("$name did not match the run or sample ID regexs");
    }
    print $run_meta_info. $/;

    return $run_meta_info;
}

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

