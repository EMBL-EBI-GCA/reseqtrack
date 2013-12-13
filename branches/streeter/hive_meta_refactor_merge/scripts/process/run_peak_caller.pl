#!/usr/bin/env perl

use strict;
use warnings;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list);
use ReseqTrack::Tools::RunPeakCall::Fseq;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use ReseqTrack::Tools::GeneralUtils
  qw(execute_system_command execute_pipe_system_command);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5);
use ReseqTrack::Tools::AttributeUtils qw(create_attribute_for_object);
use ReseqTrack::Tools::Loader::File;

use Statistics::Descriptive;

use File::Copy "move";
use File::Basename qw(fileparse);

use Data::Dumper;
use Getopt::Long;

$| = 1;

my ( $dbhost, $dbuser, $dbpass, $dbport, $dbname );
my ( $name, $type_input, $type_output, $output_dir, $directory_layout,
  $gzip_output );
my ( $store, $verbose, $program, $disable_md5 );
my $control_type;
my $strip_duplicates;
my $do_peak_stats = 1;
my $metric_type   = 'REGION_METRICS';
my $module_path   = 'ReseqTrack::Tools::RunPeakCall';
my $module_name;
my $host_name               = '1000genomes.ebi.ac.uk';
my $bam_to_bed_path         = 'bamToBed';
my $samtools_path           = 'samtools';
my $intersect_bed_path      = 'intersectBed';
my $control_experiment_type = 'ChIP-Seq Input';
my $save_temp_files         = 0;
my $run_id_regex            = '[ESD]RR\d{6}';
my $sample_id_regex         = '[ESD]RS\d{6}';
my $require_control         = 0;
my $preserve_paths          = 0;
;
my $fragment_size_stat_name;
my %options;

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
  'verbose!'                  => \$verbose,
  'program=s'                 => \$program,
  'disable_md5!'              => \$disable_md5,
  'options=s'                 => \%options,
  'bam_to_bed_path=s'         => \$bam_to_bed_path,
  'strip_duplicates!'         => \$strip_duplicates,
  'samtools_path=s'           => \$samtools_path,
  'module_path=s'             => \$module_path,
  'module_name=s'             => \$module_name,
  'intersect_bed_path=s'      => \$intersect_bed_path,
  'peak_stats!'               => \$do_peak_stats,
  'control_type=s'            => \$control_type,
  'control_experiment_type=s' => \$control_experiment_type,
  'save_temp_files=s'         => \$save_temp_files,
  'run_id_regex=s'            => \$run_id_regex,
  'sample_id_regex=s'         => \$sample_id_regex,
  'require_control!'          => \$require_control,
  'preserve_paths'            => \$preserve_paths,
  'fragment_size_stat_name=s'  => \$fragment_size_stat_name,
);

throw("Must specify an output directory ")
  if ( !$output_dir && !$preserve_paths );
throw("Must specify an output type") if ( !$type_output );

my $peak_call_module = load_module( $module_path, $module_name );
my $allowed_options = get_allowed_options($peak_call_module);
foreach my $option ( keys %options ) {
  throw( "Don't recognise option $option. Acceptable options are: "
      . join( ' ', @$allowed_options ) )
    if ( !grep { $option eq $_ } @$allowed_options );
}

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
throw("Failed to find a collection for "
    . $name . " "
    . $type_input
    . " from "
    . $dbname )
  unless ($collection);

my $input_files     = $collection->others;
my $job_name        = $name;
my @input_filepaths = map { $_->{'name'} } @$input_files;

my $run_meta_info = get_run_meta_info( $name, $rmia );
if ($preserve_paths) {
  my ( $filename, $basepath, $suffix ) = fileparse( $input_filepaths[0] );
  print Dumper( [ $filename, $basepath, $suffix ] );
  $output_dir = $basepath;
  $job_name   = $filename;
  $job_name =~ s/\.bam$//;
}
else {
  $output_dir =
    create_directory_path( $run_meta_info, $directory_layout, $output_dir )
    if ($directory_layout);
}

my $control_filepaths =
  get_control_file_paths( $name, $run_meta_info, $control_experiment_type,
  $control_type )
  if ($control_type);

throw("No control found for $name $type_input")
  if ( $require_control && (!$control_filepaths || !@$control_filepaths ));

print 'Found input files ' . join( ', ', @input_filepaths ) . $/ if $verbose;
print 'Found control files ' . join( ', ', @$control_filepaths ) . $/
  if $verbose;


print STDOUT join( "\t", "Using input:",   @input_filepaths ) . $/;
print STDOUT join( "\t", "Using control:", @$control_filepaths ) . $/ if $control_filepaths;



my %peak_caller_args = (
  -input_files              => \@input_filepaths,
  -control_files            => $control_filepaths,
  -program                  => $program,
  -working_dir              => $output_dir,
  -options                  => \%options,
  -echo_cmd_line            => $verbose,
  -job_name                 => $job_name,
  -strip_duplicates         => $strip_duplicates,
  -bam_to_bed_path          => $bam_to_bed_path,
  -samtools_path            => $samtools_path,
  -save_files_from_deletion => $save_temp_files,
);

if ($fragment_size_stat_name) {
  my $stats = $collection->statistics();
  my $fragment_size;
  for my $stat (@$stats) {
    if ($stat->attribute_name eq $fragment_size_stat_name) {
      # ppqt can find multiple peaks. the first is the most likely, so we use that one
            ($fragment_size) = split /,/, $stat->attribute_value();
    }
  }
  throw("Could not find statistic $fragment_size_stat_name for collection $name $type_input") unless $fragment_size;

  $peak_caller_args{-fragment_size} = $fragment_size;
}

print Dumper( \%peak_caller_args );

# run analysis
my $peak_caller = $peak_call_module->new(%peak_caller_args);

$peak_caller->run();
my @files = @{ $peak_caller->output_files };

for my $output (@files) {
  throw("Cannot find output file: $output") if ( !-e $output );
}

print "Peak calling complete$/" if $verbose;

my %peak_stats = gather_peak_stats( $peak_caller->bed_file, \@input_filepaths )
  if ( $do_peak_stats && $peak_caller->bed_file );

if ($gzip_output) {
  for ( my $i = 0 ; $i < scalar(@files) ; $i++ ) {
    execute_system_command("gzip -f $files[$i]");
    $files[$i] .= '.gz';
  }
}
$db->dbc->disconnect_when_inactive(0);
my $host = get_host_object( $host_name, $db );

foreach my $path (@files) {
  throw("database already has file with name $path")
    if ( $store && $fa->fetch_by_name($path) );
}

my $file_objects =
  create_objects_from_path_list( \@files, $type_output, $host );
if ( !$disable_md5 ) {
  foreach my $file_object (@$file_objects) {
    $file_object->md5( run_md5( $file_object->name ) );
  }
}
my $collection = ReseqTrack::Collection->new(
  -name   => $name,
  -type   => $type_output,
  -others => $file_objects
);

if ($do_peak_stats) {
  my @stats;

  while ( my ( $key, $value ) = each %peak_stats ) {
    push @stats, create_attribute_for_object( $collection, $key, $value )
      if ($key);
  }
  $collection->attributes( \@stats );

  if ($metric_type) {
    my $metrics_file = $output_dir . '/' . $name . '.region_metrics';
    open METRICS, '>', $metrics_file
      or throw("Could not write to $metrics_file: $!");

    my @keys = sort keys %peak_stats;
    print METRICS join( "\t", @keys ) . $/;
    for my $key (@keys) {
      print METRICS $peak_stats{$key} . "\t";
    }
    print METRICS $/;

    close METRICS;

    create_output_records( $name, $metric_type, [$metrics_file] );
  }

  $ca->store($collection) if ($store);
}

$ca->store($collection) if ($store);

sub load_module {
  my ( $module_path, $module_name ) = @_;
  my $full_module_path = join( '::', $module_path, $module_name );
  my $file = "$full_module_path.pm";
  $file =~ s{::}{/}g;
  eval { require "$file"; };
  if ($@) {
    throw("cannot load $file: $@");
  }
  return $full_module_path;
}

sub get_allowed_options {
  my $alignment_module    = shift;
  my $default_options_sub = '&' . $alignment_module . '::DEFAULT_OPTIONS';
  return [] if ( !eval "exists $default_options_sub" );
  my $default_options = eval "$default_options_sub";
  my @allowed_options = keys %$default_options;
  return \@allowed_options;
}

sub get_total_reads {
  my ($input_filepaths) = @_;

  my $read_count = 0;

  for my $f (@$input_filepaths) {
    my $cmd;
    if ( $f =~ m/\.bed$/ ) {
      $cmd = "cat $f | grep -v \\# | wc -l";
    }
    elsif ( $f =~ m/\.bam$/ ) {
      $cmd = "samtools view -c $f";
    }
    else {
      throw("Cannot get read count for $f");
    }

    print "Executing pipe command: $cmd $/" if $verbose;
    $read_count = +execute_pipe_system_command($cmd);

    print "Read count $read_count$/" if $verbose;
  }

  return $read_count;
}

sub get_reads_in_peaks {
  my ( $input_filepaths, $bed_file_path, $intersect_bed_path ) = @_;

  my $reads_in_peaks = 0;

  for my $f (@$input_filepaths) {
    my $cmd;
    if ( $f =~ m/\.bed$/ ) {
      $cmd =
"cat $bed_file_path | grep -v \\# | $intersect_bed_path -a $f -b stdin -bed | wc -l";
    }
    elsif ( $f =~ m/\.bam$/ ) {
      $cmd =
"cat $bed_file_path | grep -v \\# | $intersect_bed_path -abam $f -b stdin -bed | wc -l";
    }
    else {
      throw("Cannot get intersection count for $f");
    }
    print "Executing pipe command: $cmd $/" if $verbose;
    $reads_in_peaks = +execute_pipe_system_command($cmd);

    print "Reads in peaks $reads_in_peaks$/" if $verbose;
  }

  return $reads_in_peaks;
}

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

sub get_control_file_paths {
  my ( $name, $run_meta_info, $control_experiment_type, $control_type ) = @_;

  my @control_filepaths;

  my $ctrl_stmt = <<'QUERY_END';
select distinct c.collection_id
from run_meta_info t
join run_meta_info ctrl on t.sample_id = ctrl.sample_id
join collection c on c.name in( ctrl.run_id, ctrl.experiment_id )
where ? in( t.run_id, t.experiment_id )
and c.type = ?
and ctrl.experiment_type = 'ChIP-Seq Input'
QUERY_END

  my $query_sth = $db->dbc->prepare($ctrl_stmt);
  $query_sth->execute( $name, $control_type ) or die $query_sth->errstr;

  while ( my $rs = $query_sth->fetchrow_arrayref ) {
    my ($c_id) = @$rs;
    my $control_collection = $ca->fetch_by_dbID($c_id);
    my @file_paths = map { $_->name } @{ $control_collection->others };
    push @control_filepaths, @file_paths;
  }

  return \@control_filepaths;
}

sub gather_peak_stats {
  my ( $bed_file, $input_filepaths ) = @_;

  my $total_reads = get_total_reads($input_filepaths);
  my $reads_in_peaks =
    get_reads_in_peaks( $input_filepaths, $bed_file, $intersect_bed_path );

  if ($total_reads) {
    $peak_stats{peak_enrichment} = ( $reads_in_peaks / $total_reads );
  }
  else {
    $peak_stats{peak_enrichment} = undef;
  }

  my $stats = Statistics::Descriptive::Full->new();
  open( my $bed_fh, '<', $bed_file ) or die "Could not open $bed_file";
  while (<$bed_fh>) {
    next if (m/^#/);
    chomp;
    my ( $seq, $start, $end ) = split /\t/;
    my $length = $end - $start + 1;
    $stats->add_data($length);
  }
  close $bed_fh;

  $peak_stats{peak_count}           = $stats->count();
  $peak_stats{peak_length_mean}     = $stats->mean();
  $peak_stats{peak_length_median}   = $stats->median();
  $peak_stats{peak_length_std_dev}  = $stats->standard_deviation();
  $peak_stats{peak_length_variance} = $stats->variance();

  map { print "Stat: $_ " . $peak_stats{$_} . $/ } keys %peak_stats
    if ($verbose);

  return %peak_stats;
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

