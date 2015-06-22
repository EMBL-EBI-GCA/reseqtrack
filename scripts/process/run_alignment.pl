#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::RunMetaInfoUtils qw(create_directory_path);
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $name;
my $type;
my $output_dir;
my $reference;
my $ref_samtools_index;
my $module_path = 'ReseqTrack::Tools::RunAlignment';
my $module_name;
my $module_constructor_args;
my $program_file;
my $samtools;
my $host_name = '1000genomes.ebi.ac.uk';
my $store;
my $update;
my $bam_type = 'BAM';
my $bai_type = 'BAI';
my $merge_bams = 1;
my $sort_bams = 1;
my $index_bams = 1;
my $convert_sam_to_bam = 1;
my $directory_layout;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'name=s' => \$name,
  'type=s' => \$type,
  'bam_type=s' => \$bam_type,
  'bai_type=s' => \$bai_type,
  'output_dir=s' => \$output_dir,
  'reference=s' => \$reference,
  'ref_samtools_index=s' => \$ref_samtools_index,
  'module_path=s' => \$module_path,
  'module_name=s' => \$module_name,
  'module_constructor_args=s' => \$module_constructor_args,
  'program_file=s' => \$program_file,
  'samtools=s' => \$samtools,
  'host_name=s' => \$host_name,
  'store!' => \$store,
  'update!' => \$update,
  'merge_bams=s' => \$merge_bams,
  'index_bams=s' => \$index_bams,
  'sort_bams=s' => \$sort_bams,
  'convert_sam_to_bam=s' => \$convert_sam_to_bam,
  'directory_layout=s' => \$directory_layout,
    );

$directory_layout ||= 'population/sample_name/' . $module_name;

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );
$db->dbc->disconnect_when_inactive(1);
my $ca = $db->get_CollectionAdaptor;


my $collection = $ca->fetch_by_name_and_type($name, $type);

throw("Failed to find a collection for ".$name." ".$type." from ".$dbname) 
    unless($collection);

$name =~ /[ES]RR\d+/ or throw("could not get run_id from $name");
my $run_id = $&;
my $rmi_a = $db->get_RunMetaInfoAdaptor;
my $run_meta_info = $rmi_a->fetch_by_run_id($run_id);
throw("Failed to find run_meta_info for $run_id from $dbname")
    if(!$run_meta_info);
my $paired_length = $run_meta_info->paired_length;
my $full_output_dir = create_directory_path($run_meta_info, $directory_layout, $output_dir);

my $input_files = $collection->others;
my @input_file_paths = map {$_->{'name'}} @$input_files;

my $alignment_module = $module_path."::".$module_name;
my $constructor_hash = parameters_hash($module_constructor_args);
$constructor_hash->{-input_files} = \@input_file_paths;
$constructor_hash->{-program} = $program_file;
$constructor_hash->{-working_dir} = $full_output_dir;
$constructor_hash->{-reference} = $reference;
$constructor_hash->{-ref_samtools_index} = $ref_samtools_index;
$constructor_hash->{-samtools} = $samtools;
$constructor_hash->{-job_name} = $name;
$constructor_hash->{-paired_length} = $paired_length;
$constructor_hash->{-merge_bams} = $merge_bams;
$constructor_hash->{-sort_bams} = $sort_bams;
$constructor_hash->{-index_bams} = $index_bams;
$constructor_hash->{-convert_sam_to_bam} = $convert_sam_to_bam;

my $run_alignment = setup_alignment_module($alignment_module, $constructor_hash);

$run_alignment->run;

my $bam_filepaths = $run_alignment->output_bam_files;
my $bai_filepaths = $run_alignment->output_bai_files;

if($store){
  my $host = get_host_object($host_name, $db);
  my $bam_files = create_objects_from_path_list($bam_filepaths, $bam_type, $host);
  my $bai_files = create_objects_from_path_list($bai_filepaths, $bai_type, $host);

  my $fa = $db->get_FileAdaptor;
  my $ca = $db->get_CollectionAdaptor;
  foreach my $file(@$bam_files, @$bai_files){
    my $md5 = run_md5($file->name);
    $file->md5($md5);
    $fa->store($file, $update);
  }
  my $collection = ReseqTrack::Collection->new(
      -name => $name, -type => $bam_type,
      -table_name => 'file', -others => $bam_files);
  $ca->store($collection);
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
