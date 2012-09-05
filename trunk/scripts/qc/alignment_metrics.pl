#!/sw/arch/bin/perl -w

use strict;
use warnings;

use ReseqTrack::Tools::RunPicard;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::StatisticsUtils;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::DBSQL::DBAdaptor;

use Getopt::Long;
use Data::Dumper;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;

my $java_path = '/usr/bin/java';
my $picard_dir;
my $file_name;
my $jvm_options = '-Xmx2g -Xms2g';
my $reference_sequence;
my $host_name = '1000genomes.ebi.ac.uk';
my $keep_metrics_file = 0;
my $metrics_file_type = 'BAM_METRICS';

&GetOptions( 
	'dbhost=s' => \$dbhost,
	'dbname=s' => \$dbname,
	'dbuser=s' => \$dbuser,
	'dbpass=s' => \$dbpass,
	'dbport=s' => \$dbport,
	'host=s' => \$host_name,
	'java_path=s' => \$java_path,
	'picard_dir=s' => \$picard_dir,
	'jvm_options=s' => \$jvm_options,
	'reference_sequence=s' => \$reference_sequence,
	'keep_metrics_file!' => \$keep_metrics_file,
	'file=s' => \$file_name,
	'metrics_file_type=s' => \$metrics_file_type,
);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );
$db->dbc->disconnect_when_inactive(1);
my $host = get_host_object($host_name,$db);

my $file_adaptor = $db->get_FileAdaptor;
my $file = $file_adaptor->fetch_by_name($file_name);

my $options = {
	"VALIDATION_STRINGENCY" => "LENIENT",
	"REFERENCE_SEQUENCE" => $reference_sequence,
};

my %picard_config = (
	-program		=> $java_path,
	-picard_dir		=> $picard_dir,
	-jvm_options	=> $jvm_options,
	-options		=> $options,
	-input_files	=> [$file_name],
);

my $picard = ReseqTrack::Tools::RunPicard->new(%picard_config);
my $alignment_metrics = $picard->run_alignment_metrics;
my $metrics_file = $picard->output_files->[0];

my $statistics = $file->statistics;

for my $metrics (@$alignment_metrics) {
	my $category = $metrics->{'CATEGORY'};
	while (my ($key, $value) = each %$metrics) {
		next if $key eq 'CATEGORY';
		push @$statistics, create_statistic_for_object($file,$category.'_'.$key,$value) if (defined $value);
	}	 
}

$file->uniquify_statistics($statistics);
$file_adaptor->store_statistics($file);

if ($keep_metrics_file){
	my $metrics_file = create_object_from_path($metrics_file,$metrics_file_type,$host);
	
	my $md5 = run_md5($metrics_file->full_path);
	$metrics_file->md5($md5);

	my $size = -s $metrics_file->full_path;
	$metrics_file->size($size);
	
	$file_adaptor->store($metrics_file);
}
else {
	$picard->files_to_delete($metrics_file);
}
$picard->delete_files;