#!/usr/bin/env perl
use strict;
use Getopt::Long;
use Data::Dumper;
use warnings;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use Statistics::Descriptive;
use ReseqTrack::Tools::GeneralUtils
  qw(execute_system_command execute_pipe_system_command);
use autodie;
use ReseqTrack::Tools::AttributeUtils qw(create_attribute_for_object);
use ReseqTrack::Tools::Loader::File;
use File::Basename qw(fileparse);

$ENV{"LD_LIBRARY_PATH"} = '';

my ( $dbhost, $dbuser, $dbpass, $dbport, $dbname );
my ( $name, $peak_type, $bam_type, $metrics_type, );
my ( $samtools_path, $bedtools_path );
my ( $print_out, $store_file, $store_attributes, $clobber );
my $host_name = '1000genomes.ebi.ac.uk';
my $peaks_re  = '\.\d{8}\.bed\.gz$';

&GetOptions(
    'dbhost=s'          => \$dbhost,
    'dbname=s'          => \$dbname,
    'dbuser=s'          => \$dbuser,
    'dbpass=s'          => \$dbpass,
    'dbport=s'          => \$dbport,
    'name=s'            => \$name,
    'peak_type=s'       => \$peak_type,
    'bam_type=s'        => \$bam_type,
    'metrics_type=s'    => \$metrics_type,
    'samtools_path=s'   => \$samtools_path,
    'bedtools_path=s'   => \$bedtools_path,
    'print_out!'        => \$print_out,
    'store_file!'       => \$store_file,
    'store_attributes!' => \$store_attributes,
    'clobber!'          => \$clobber,
    'host_name=s'       => \$host_name,
    'peaks_re=s'        => \$peaks_re,
);

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
my $ca = $db->get_CollectionAdaptor;

my $peaks_file = get_file( $ca, $name, $peak_type, $peaks_re, 1, );
my $bam_file   = get_file( $ca, $name, $bam_type,  '\.bam$',  1, );
my $metrics_file_name =
  metrics_file_name( $peaks_file, $clobber, $metrics_type )
  if ($store_file);

my %peak_stats =
  gather_peak_stats( $peaks_file, $bam_file, $samtools_path, $bedtools_path );

if ($print_out) { print Dumper( \%peak_stats ); }

if ($store_file) {

    open METRICS, '>', $metrics_file_name;
    my @keys = sort keys %peak_stats;
    print METRICS join( "\t", @keys ) . $/;
    for my $key (@keys) {
        print METRICS $peak_stats{$key} . "\t";
    }
    print METRICS $/;
    close METRICS;

    #create a File loader object for the files
    my $loader = ReseqTrack::Tools::Loader::File->new(
        -file            => [$metrics_file_name],
        -do_md5          => 1,
        -hostname        => $host_name,
        -db              => $db,
        -assign_types    => 0,
        -check_types     => 0,
        -type            => $metrics_type,
        -update_existing => 1,
    );

    #Load the files into the database
    $loader->process_input();
    $loader->create_objects();
    $loader->sanity_check_objects();
    my $file_objects = $loader->load_objects();

    my $metrics_coll = $ca->fetch_by_name_and_type( $name, $metrics_type );
    if ($metrics_coll) {
        $metrics_coll->others($file_objects);
    }
    else {
        $metrics_coll = ReseqTrack::Collection->new(
            -name       => $name,
            -type       => $metrics_type,
            -others     => $file_objects,
            -table_name => 'file',
        );

    }
    $ca->store($metrics_coll);
}

if ($store_attributes) {
    my $peak_collection = $ca->fetch_by_name_and_type( $name, $peak_type );

    my @stats;

    while ( my ( $key, $value ) = each %peak_stats ) {
        push @stats,
          create_attribute_for_object( $peak_collection, $key, $value )
          if ( defined $key && defined $value );
    }
    $peak_collection->attributes( \@stats );
    $ca->store($peak_collection);
}

sub metrics_file_name {
    my ( $peaks_file, $clobber,  $metrics_type ) = @_;
    my ( $filename,   $basepath, $suffix )       = fileparse($peaks_file);
    my $metrics_file_name = $basepath . '/' . $name . '.' . lc($metrics_type);
    die("Metrics file already exists: $metrics_file_name")
      if ( -e $metrics_file_name && !$clobber );
    return $metrics_file_name;
}

sub get_file {
    my ( $ca, $name, $type, $grep_block, $required ) = @_;
    my $coll = $ca->fetch_by_name_and_type( $name, $type );

    die "No collection found for $name $type" unless $coll;

    my @files = grep /$grep_block/, map { $_->name } @{ $coll->others };

    die "No file found for $name $type $grep_block" if ( $required && !@files );
    die "Multiple files found for $name $type $grep_block" if ( @files > 1 );

    return pop @files;
}

sub gather_peak_stats {
    my ( $bed_file, $bam_file, $samtools_path, $bedtools_path ) = @_;

    my %peak_stats;

    $peak_stats{total_reads} = get_total_reads( $bam_file, $samtools_path );
    $peak_stats{reads_in_peaks} =
      get_reads_in_peaks( $bam_file, $bed_file, $bedtools_path );

    if ( $peak_stats{total_reads} ) {
        $peak_stats{peak_enrichment} =
          ( $peak_stats{reads_in_peaks} / $peak_stats{total_reads} );
    }
    else {
        $peak_stats{peak_enrichment} = undef;
    }

    my $stats = Statistics::Descriptive::Full->new();

    my $bed_fh;
    if ( $bed_file =~ m/\.gz$/ ) {
        open $bed_fh, '-|', 'gzip', '-dc', $bed_file;
    }
    else {
        open( $bed_fh, '<', $bed_file );
    }
    while (<$bed_fh>) {
        next if (m/^#/);
        chomp;
        my ( $seq, $start, $end ) = split /\t/;
        my $length = $end - $start + 1;
        $stats->add_data($length);
    }
    close $bed_fh;

    $peak_stats{region_count}           = $stats->count();
    $peak_stats{region_length_mean}     = $stats->mean();
    $peak_stats{region_length_median}   = $stats->median();
    $peak_stats{region_length_std_dev}  = $stats->standard_deviation();
    $peak_stats{region_length_variance} = $stats->variance();
    $peak_stats{region_length_total}    = $stats->sum();

    return %peak_stats;
}

sub get_total_reads {
    my ( $bam_file, $samtools_path ) = @_;

    my $read_count =
      execute_pipe_system_command("$samtools_path view -c $bam_file");

    return $read_count;
}

sub get_reads_in_peaks {
    my ( $bam_file, $bed_file, $bedtools_path ) = @_;

    my $reads_in_peaks = execute_pipe_system_command(
        "$bedtools_path intersect -a $bam_file -b $bed_file -bed | wc -l");

    return $reads_in_peaks;
}
