#!/sw/arch/bin/perl -w

use strict;
use Getopt::Long;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::SequenceIndexUtils;

$| = 1;

my $dbhost;
my $dbname;
my $dbuser;
my $dbport;
my $dbpass;
my $run_id; 
my $type = 'ARCHIVE_FASTQ';

&GetOptions( 
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'run_id=s' => \$run_id,
  'type=s' => \$type,
    );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $rmi_a = $db->get_RunMetaInfoAdaptor;
my $meta_info = $rmi_a->fetch_by_run_id($run_id);


my $ca = $db->get_CollectionAdaptor;

my $collection = $ca->fetch_by_name_and_type($run_id, $type);

if(!$collection){
  throw("Failed to fetch a collection for ".$run_id." ".$type." from ".$dbname);
}


my $files = $collection->others;
throw("Collection ".$collection->name." ".$collection->type." in ".$dbname." has no files") 
    if(!$files || @$files == 0);
my @paths;
foreach my $file(@$files){
  push(@paths, $file->name);
}

my @problems;
foreach my $path(@paths){
  push(@problems, $path." doesn't exist") unless(-e $path);
}

my ($mate1, $mate2, $frag) = assign_files(\@paths);


if($meta_info->library_layout eq 'SINGLE'){
  if($mate1 || $mate2){
    my $problem = "Have mate pair files for ".$meta_info->run_id." ".$meta_info->library_layout;
    $problem .= " ".$mate1." " if($mate1);
    $problem .= " and " if($mate1 && $mate2);
    $problem .= " ".$mate2." " if($mate2);
    $problem .= " are defined";
    push(@problems, $problem);
  }
  if(!$frag){
    my $problem .= "Don't have a fragment/single ended fastq for ".$meta_info->run_id." ".$meta_info->library_layout;
    push(@problems, $problem);
  }
}elsif($meta_info->library_layout eq 'PAIRED'){
  if(!$mate1 || !$mate2){
    my $problem = $meta_info->run_id." is missing mate pair files ";
    $problem .= "MATE 1 ".$mate1." " if($mate1);
    $problem .= "MATE 2 ".$mate2." " if($mate2);
  }
}else{
  push(@problems, "Don't know how to treat ".$meta_info->library_layout);
}

print "Here are the problems\n";
foreach my $problem(@problems){
  print STDERR $problem."\n";
}

if(@problems >=1){
  throw("There are ".@problems." with ".$run_id);
}
