#!/sw/arch/bin/perl -w

use strict;

use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::RunMetaInfoUtils;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;

$| = 1;

my $index_file;
my $output_dir;
my $output_name;
my $help;
my $clobber;
my $collection_type = 'STUDY_TYPE';
my $collection_name;
my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;

&GetOptions( 
	    'dbhost=s'       => \$dbhost,
	    'dbname=s'       => \$dbname,
	    'dbuser=s'       => \$dbuser,
	    'dbpass=s'       => \$dbpass,
	    'dbport=s'       => \$dbport,
	    'index_file=s' => \$index_file,
	    'output_dir=s' => \$output_dir,
	    'output_name=s' => \$output_name,
	    'help!' => \$help,
	    'collection_type:s' => \$collection_type,
	    'collection_name:s' => \$collection_name,
	    'clobber!' => \$clobber,
    );

if($help){
  useage();
}

throw("Can't dump stats without an index file you must pass on in using ".
      "-index_file and this file must exist")
    unless($index_file && -e $index_file);

throw("Can't output to ".$output_dir." directory it doesn't exist")
    unless($output_dir && -d $output_dir);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );


unless($output_name){
  my $index_name = basename($index_file);
  $output_name = $index_name;
  $output_name =~ s/sequence\.index/sequence\.index\.stats/;
}

my $full_output_path = $output_dir."/".$output_name;
if(-e $full_output_path){
  throw("Can't use ".$full_output_path." as it already exists") unless($clobber);
}

my @rmis;


if($collection_name && $collection_type){
  my $ca = $db->get_CollectionAdaptor;
  my $collection =  $ca->fetch_by_name_and_type($collection_name, $collection_type);
  @rmis = @{$collection->others};
}else{
  my $rmia = $db->get_RunMetaInfoAdaptor;
  @rmis = @{$rmia->fetch_all};
}

my %run_id_hash;

foreach my $rmi(@rmis){
  $run_id_hash{$rmi->run_id} = 1;
  #print $rmi->run_id . "\t" . $rmi->study_id . "\t" . $collection_name . "\n";  
}

my $population_rules = $db->get_PopulationRuleAdaptor->fetch_all_in_order();

my $withdrawn_hash = get_withdrawn_summary($index_file, \%run_id_hash);
my $study_hash = get_sequence_index_stats($index_file, 3, 24, \%run_id_hash, $population_rules);
my $study_desc = get_study_descriptions($index_file);
my $center_hash = get_index_group_stats($index_file, 5, 3, 24, \%run_id_hash, $population_rules);
my $pop_hash = get_sequence_index_stats($index_file, 10, 24, \%run_id_hash, $population_rules);
my $sample_hash = get_sequence_index_stats($index_file, 9, 24, \%run_id_hash, $population_rules);

open(FH, ">".$full_output_path) or throw("Failed to open ".$full_output_path." $!");
#print withdrawn summary
my $maxreason = 0; 
my $maxnumber = 0;
foreach my $key(keys%$withdrawn_hash){
  my $number = $withdrawn_hash->{$key};
  $maxreason = length($key) if(length($key) > $maxreason);
  $maxnumber = length($number) if(length($number) > $maxnumber);
}
print FH "WITHDRAWN SUMMARY\n\n";
my $string = sprintf  ("%-${maxreason}s\t%-${maxnumber}s\n","REASON","NUMBER");
print FH $string;
my @keys = sort{$a cmp $b} keys(%$withdrawn_hash);
foreach my $key(@keys){
  my $number = $withdrawn_hash->{$key};
  $string = sprintf ("%-${maxreason}s\t%-${maxnumber}s\n", $key, $number);
  print FH $string;
}
print FH "\n";
my $maxstudy = 0;
my $maxdesc = 0;
my $maxbp = 0;
my $maxcoverage = 0;
@keys = sort{$a cmp $b} keys(%$study_hash);
foreach my $key(@keys){
  my $desc = $study_desc->{$key};
  my $count = $study_hash->{$key};
  my $coverage = calculate_coverage($count);
  $maxstudy = length($key) if(length($key) > $maxstudy);
  $maxdesc = length($desc) if(length($desc) > $maxdesc);
  $maxbp = length($count) if(length($count) > $maxbp);
  $maxcoverage = length($coverage) if(length($coverage) > $maxcoverage);
}
print FH "STUDY SUMMARY\n\n";
$string = sprintf("%-${maxstudy}s\t%-${maxdesc}s\t%-${maxbp}s\t%-${maxcoverage}s\n",
            "STUDY ID", "STUDY NAME", "BASE COUNT", "COVERAGE");
print FH $string;
foreach my $key(@keys){
  my $desc = $study_desc->{$key};
  my $count = $study_hash->{$key};
  my $coverage = calculate_coverage($count);
  $string = sprintf("%-${maxstudy}s\t%-${maxdesc}s\t%-${maxbp}s\t%-${maxcoverage}s\n",
                    $key, $desc, $count, $coverage);
  print FH $string;
}
print FH "\n";
my $maxpop = 0;
foreach my $key(keys(%$pop_hash)){
  my $count = $pop_hash->{$key};
  my $coverage = calculate_coverage($count);
  $maxpop = length($key) if(length($key) > $maxpop);
  $maxbp = length($count) if(length($count) > $maxbp);
  $maxcoverage = length($coverage) if(length($coverage) > $maxcoverage);
}
print FH "POPULATION SUMMARY\n\n";
$string = sprintf  ("%-${maxpop}s\t%-${maxbp}s\t%-${maxcoverage}s\n",
                    "POP", "BASE COUNT", "COVERAGE");
print FH $string;
@keys = sort{$a cmp $b} keys(%$pop_hash);
foreach my $key(@keys){
  my $count = $pop_hash->{$key};
  my $coverage = calculate_coverage($count);
  $string = sprintf  ("%-${maxpop}s\t%-${maxbp}s\t%-${maxcoverage}s\n",
            $key, $count, $coverage);
  print FH $string;
}
print FH "\n";
my $maxcenter = 0;
foreach my $key(keys(%$center_hash)){
  my $hash = $center_hash->{$key};
  foreach my $study(keys(%$hash)){
    my $count = $hash->{$study};
    my $coverage = calculate_coverage($count);
    $maxcenter = length($key) if(length($key) > $maxcenter);
    $maxbp = length($count) if(length($count) > $maxbp);
    $maxcoverage = length($coverage) if(length($coverage) > $maxcoverage);
  }
}
print FH "CENTER SUMMARY\n\n";
$string = sprintf("%-${maxcenter}s\t%-${maxstudy}s\t%-${maxbp}s\t%-${maxcoverage}s\n",
                  "CENTER", "STUDY", "BASE COUNT", "COVERAGE");
print FH $string;
@keys = sort{$a cmp $b} keys(%$center_hash);
foreach my $key(@keys){
  my $hash = $center_hash->{$key};
  my @sorted = sort{$a cmp $b} keys(%$hash);
  foreach my $study(@sorted){
    my $count = $hash->{$study};
    my $coverage = calculate_coverage($count);
    $string = sprintf("%-${maxcenter}s\t%-${maxstudy}s\t%-${maxbp}s\t%-${maxcoverage}s\n",
                      $key, $study, $count, $coverage);
    print FH $string;
  }
}
print FH "\n";
my $maxsamp = 0;
foreach my $key(keys(%$sample_hash)){
  my $count = $sample_hash->{$key};
  my $coverage = calculate_coverage($count);
  $maxsamp = length($key) if(length($key) > $maxsamp);
  $maxbp = length($count) if(length($count) > $maxbp);
  $maxcoverage = length($coverage) if(length($coverage) > $maxcoverage);
}
print  FH "SAMPLE SUMMARY\n\n";
$string = sprintf  ("%-${maxsamp}s\t%-${maxbp}s\t%-${maxcoverage}s\n",
            "SAMPLE", "BASE COUNT", "COVERAGE");
print FH $string;
@keys = sort{$a cmp $b} keys(%$sample_hash);
foreach my $key(@keys){
  my $count = $sample_hash->{$key};
  my $coverage = calculate_coverage($count);
  $string = sprintf  ("%-${maxsamp}s\t%-${maxbp}s\t%-${maxcoverage}s\n",
            $key, $count, $coverage);
  print FH $string;
}


=pod

=head1 NAME

ReseqTrack/scripts/run_meta_info/dump_sequence_index_stats.pl

=head1 SYNOPSIS

This script produces summary statistics based on the base count column and 
withdrawn/withdrawn comment columns 

=head1 OPTIONS

-index_file, path to index file

-output_dir, path to output directory

-output_name, name of output file, if this isn't defined the input filename is taken and
s/sequence\.index/sequence\.index\.stats/. If the output filename isn't specified and the 
input name doesn't match this pattern there may be problems

-clobber, over write an existing file

-help, print out perldocs


=head1 Examples

perl ReseqTrack/scripts/run_meta_info/dump_sequence_index_stats.pl -index /nfs/1000g-archive/vol1/ftp/sequence_indices/20091216.sequence.index 
 -output_dir /nfs/1000g-work/G1K/archive_staging/ftp/sequence_indices/

=head1 Other useful scripts

ReseqTrack/scripts/dump_sequence_index.pl

=cut

