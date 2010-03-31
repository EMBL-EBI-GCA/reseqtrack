#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileSystemUtils;
use File::Basename;
use Getopt::Long;
use ReseqTrack::Tools::RunMetaInfoUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::GeneralUtils;

$| = 1;

my $index_file;
my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $output_file;
my $help;

&GetOptions( 
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'index_file=s' => \$index_file,
  'output_file=s' => \$output_file,
  'help!' => \$help,
    );

if($help){
  useage();
}

my $lines = get_lines_from_file($index_file);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $rmia = $db->get_RunMetaInfoAdaptor;

my $rmis = $rmia->fetch_all;

my $fh;
if($output_file){
  open(FH, ">".$output_file) or throw("Failed to open ".$output_file." $!");
  $fh = \*FH;
}else{
  $fh = \*STDOUT;
}

my %sra_hash;
my $sra_total;
foreach my $rmi(@$rmis){
  next unless($rmi->status eq 'public');
  next if($rmi->study_id eq 'SRP000032' || $rmi->study_id eq 'SRP000033');
  $sra_hash{$rmi->center_name} += $rmi->archive_base_count;
  $sra_total += $rmi->archive_base_count;
}

print $fh "SRA Statistics in Gb\n";
my @center_names = sort{$a cmp $b} keys(%sra_hash);
foreach my $center_name(@center_names){
  print $fh $center_name.", ".calculate_gigabase($sra_hash{$center_name})."\n";
}
print $fh "Total SRA Base count, ".calculate_gigabase($sra_total)."\n"; 

my %index_withdrawn_accessions;
my %index_not_available_accessions;
my %index_accessions;
my %index_samples;
my %index_pop_sample_count;
my %index_pop_accession_count;
my %index_sample_base_count;
my %index_instrument_base_count;
my %index_population_base_count;

my $total_index_base_count;
my $index_name = basename($index_file);
$index_name =~ /^(\d+)\.sequence\.index/;
my $index_date = $1;
foreach my $line(@$lines){
  next if($line =~ /SUBMISSION/);
  my @values = split /\t/, $line;
  if($values[20]){
    if($values[22] eq 'NOT YET AVAILABLE FROM ARCHIVE'){
      $index_not_available_accessions{$values[2]}++;
    }else{
      $index_withdrawn_accessions{$values[2]}++;
    }
    next;
  }
  unless($values[25] eq 'low coverage'){
    next;
  }
  if($values[24] eq 'not available'){
    throw("Don't have stats for ".$line);
  }
  $index_accessions{$values[2]} = 1;
  $index_samples{$values[9]}{$values[2]} = 1;
  $index_sample_base_count{$values[9]} += $values[24];
  $index_instrument_base_count{$values[12]} += $values[24];
  $index_population_base_count{$values[10]} += $values[24];
  $index_pop_sample_count{$values[10]}{$values[9]} = 1;
  $index_pop_accession_count{$values[10]}{$values[2]} = 1;
  $total_index_base_count += $values[24];
                         
}

my $index_samples_over_12;
foreach my $index_sample(keys(%index_sample_base_count)){
  my $gigabase = calculate_gigabase($index_sample_base_count{$index_sample}, 1);
  $index_samples_over_12++ if($gigabase >= 12);
}

print $fh "\nDCC index summary\n";
print $fh "Last Sequence Index, ".$index_date."\n";
print $fh "\t# of accessions, ".keys(%index_accessions)."\n";
print $fh "\t# of samples, ".keys(%index_samples)."\n";
print $fh "\t# of samples > 4x, ".$index_samples_over_12."\n";
print $fh "\tTotal Gb in index, ".calculate_gigabase($total_index_base_count)."\n";
print $fh "\n";
print $fh "Population Summary, accession count, sample count, base count in GB\n";
my @sorted_pops = sort{$a cmp $b} keys(%index_population_base_count);
foreach my $pop(@sorted_pops){
  print $fh $pop.", ".keys(%{$index_pop_accession_count{$pop}}).", ".keys(%{$index_pop_sample_count{$pop}}).", ".calculate_gigabase($index_population_base_count{$pop})."\n";
}
print $fh "\n";
print $fh "Platform Summary in Gb\n";
my @sorted_instruments = sort{$a cmp $b} keys(%index_instrument_base_count);
foreach my $instrument(@sorted_instruments){
  print $fh $instrument.", ".calculate_gigabase($index_instrument_base_count{$instrument})."\n";
}


=pod

=head1 NAME

ReseqTrack/scripts/run_meta_info/dcc_summary_stats.pl

=head1 SYNOPSIS

This script uses an index file and the archive base counts in the run meta info
table to give a summary of sequence available in the system. The stats are output
in csv format.

=head1 OPTIONS

-index_file, path to index file

-output_file, path to output file, if no file is specified the script prints to
STDOUT

-dbhost, the name of the mysql-host

-dbname, the name of the mysql database

-dbuser, the name of the mysql user

-dbpass, the database password if appropriate

-dbport, the port the mysql instance is running on, this defaults to 4197 the 
    standard port for mysql-g1kdcc.ebi.ac.uk

-help, print out perldocs


=head1 Examples

perl ReseqTrack/scripts/run_meta_info/dcc_summary_stats.pl -dbhost myhost -dbuser myuser -dbport 3306 -dbname mydatabase -index 20100311.sequence.index -output_file /path/to/output_file
=head1 Other useful scripts

ReseqTrack/scripts/dump_sequence_index.pl

=cut
