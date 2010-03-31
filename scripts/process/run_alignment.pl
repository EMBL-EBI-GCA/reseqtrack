#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
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
my $reference;
my $module_path = "ReseqTrack::Tools::BatchSubmission";
my $module_name;
my $module_constructor_args;
my $program_file;
my $samtools;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $update;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'run_id=s' => \$run_id,
  'type=s' => \$type,
  'output_dir=s' => \$output_dir,
  'reference=s' => \$reference,
  'module_path=s' => \$module_path,
  'module_name=s' => \$module_name,
  'module_constructor_args=s' => \$module_constructor_args,
  'program_file=s' => \$program_file,
  'samtools=s' => \$samtools,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'update!' => \$update,
    );

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

my $alignment_module = $module_path."::".$module_name;
my $constructor_hash = parameters_hash($module_constructor_args);
$constructor_hash->{reference} = $reference;
$constructor_hash->{input} = $collection;
$constructor_hash->{program} = $program_file;
$constructor_hash->{samtools} = $samtools;
$constructor_hash->{working_dir} = $output_dir;
$constructor_hash->{name} = $run_id;

my $run_alignment = setup_alignment_module($alignment_module, $constructor_hash);

$run_alignment->run;

my $output_bams = $run_alignment->output_files;

my $host = get_host_object($host_name, $db);
my $files = create_objects_from_path_list($output_bams, 'BAM', $host);

if($store){
  my $fa = $db->get_FileAdaptor;
  foreach my $file(@$files){
  my $md5 = run_md5($file->name);
  $file->md5($md5);
  $fa->store($file, $update);
  }
}

sub parameters_hash{
  my ($string) = @_;

  my %parameters_hash;

  if ($string) {
    if($string =~  /,/ || $string =~ /=>/){
      my @pairs = split (/,/, $string);
      foreach my $pair(@pairs){
        my ($key, $value) = split (/=>/, $pair);
        if ($key && ($value || $value == 0)) {
          $key   =~ s/^\s+//g;
          $key   =~ s/\s+$//g;
          $value =~ s/^\s+//g;
          $value =~ s/\s+$//g;
          $parameters_hash{$key} = $value;
        } else {
          $parameters_hash{$key} = 1;
        }
      }
    }else{
      $parameters_hash{'-options'} = $string;
    }
  }
  return \%parameters_hash;
}

sub setup_alignment_module{
  my ($alignment_module, $args) = @_;
  my $file = "$alignment_module.pm";
  $file =~ s{::}{/}g;
  eval {
    require "$file";
  };
  if($@){
    throw("ReseqTrack::Tools::EventPipeline::setup_batch_submission_system ".
          "Can't find $file [$@]");
  }
  my %constructor_args = %$args;
  my $object = $alignment_module->new(%constructor_args);
  return $object;
}
