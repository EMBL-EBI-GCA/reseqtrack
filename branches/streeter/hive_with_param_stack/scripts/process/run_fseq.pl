#!/usr/bin/env perl

use strict;
use warnings;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list);
use ReseqTrack::Tools::RunPeakCall::Fseq;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use ReseqTrack::Tools::GeneralUtils qw(execute_system_command);
use ReseqTrack::Tools::FileSystemUtils qw(run_md5);

use File::Copy "move";
use File::Basename qw(fileparse);

use Getopt::Long;

$| = 1;


my ($dbhost, $dbuser, $dbpass, $dbport, $dbname);
my ($name, $type_input, $type_output,$output_dir,$directory_layout,$gzip_output);
my ($store,$verbose,$program,$disable_md5);
my ($bam_to_bed_path,$samtools_path,$strip_duplicates);
my $host_name = '1000genomes.ebi.ac.uk';
my %options;
my $run_id_regex = '[ESD]RR\d{6}';
my $sample_id_regex = '[ESD]RS\d{6}';

&GetOptions( 
	'dbhost=s'      => \$dbhost,
	'dbname=s'      => \$dbname,
	'dbuser=s'      => \$dbuser,
	'dbpass=s'      => \$dbpass,
	'dbport=s'      => \$dbport,
	'name=s' => \$name,
	'type_input=s' => \$type_input,
	'type_output=s' => \$type_output,
	'output_dir=s' => \$output_dir,
	'directory_layout=s' => \$directory_layout,
	'gzip_output!' => \$gzip_output,		
	'store!' => \$store,
	'verbose!' => \$verbose,
	'program=s' => \$program,
	'disable_md5!' => \$disable_md5,
	'options=s' => \%options,
	'bam_to_bed_path=s' => \$bam_to_bed_path,
	'strip_duplicates!' => \$strip_duplicates,
	'samtools_path=s' => \$samtools_path,	
	'run_id_regex=s' => \$run_id_regex,
	'sample_id_regex=s' => \$sample_id_regex,
);

throw("Must specify an output directory") if (!$output_dir);
throw("Must specify an output type") if (!$type_output);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
);

throw("Could not create db adaptor") if (! $db);

$db->dbc->disconnect_when_inactive(1);
my $ca = $db->get_CollectionAdaptor;
my $fa = $db->get_FileAdaptor;
my $rmia = $db->get_RunMetaInfoAdaptor;



my $collection = $ca->fetch_by_name_and_type($name, $type_input);
throw("Failed to find a collection for ".$name." ".$type_input." from ".$dbname) 
    unless($collection);

my $input_files = $collection->others;
my @input_filepaths = map {$_->{'name'}} @$input_files;

print 'Found input files '.join( ', ', @input_filepaths).$/ if $verbose;

if ($directory_layout) {
  my $run_meta_info;
  if ($name =~ /$run_id_regex/) {
    $run_meta_info = $rmia->fetch_by_run_id($&);
  }
  elsif ($name =~ /$sample_id_regex/) {
    my $rmi_list = $rmia->fetch_by_sample_id($&);
    $run_meta_info = $rmi_list->[0] if (@$rmi_list);
  }

  if ($run_meta_info) {
    $output_dir = create_directory_path($run_meta_info, $directory_layout, $output_dir);
  }
}

my $fseq = ReseqTrack::Tools::RunPeakCall::Fseq->new(
	-input_files => \@input_filepaths,
	-program => $program,
	-working_dir => $output_dir,
	-options => \%options,
	-echo_cmd_line => $verbose,
	-job_name => $name,
	-strip_duplicates => $strip_duplicates,
	-bam_to_bed_path => $bam_to_bed_path,
	-samtools_path => $samtools_path,	
);

$fseq->run();

my @files = @{$fseq->output_files};

for (my $i = 0; $i < scalar(@files); $i++) {
	my ($file_name,$path) = fileparse($files[$i]);
	my $new_name = $path.$name.'_'.$file_name;
	move($files[$i], $new_name);
	
	if ($gzip_output){
		execute_system_command("gzip $new_name");
		$new_name .= '.gz';
	}
	
	$files[$i] = $new_name;
}

$db->dbc->disconnect_when_inactive(0);
my $host = get_host_object($host_name, $db);

foreach my $path (@files) {
	throw("database already has file with name $path")if ($fa->fetch_by_name($path));
}

my $file_objects = create_objects_from_path_list(\@files, $type_output, $host);
if (! $disable_md5) {
	foreach my $file_object (@$file_objects) {
		$file_object->md5( run_md5($file_object->name) );
	}
}
my $collection = ReseqTrack::Collection->new(
	-name => $name,
	-type => $type_output,
	-others => $file_objects
);

$ca->store($collection) if($store);

