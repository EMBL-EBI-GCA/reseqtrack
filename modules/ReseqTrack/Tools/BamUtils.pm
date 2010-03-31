package ReseqTrack::Tools::BamUtils;
use strict;
use warnings;

use Exporter;
use ReseqTrack::File;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use File::Copy;
use File::Basename;
use File::Find ();

use vars qw (@ISA  @EXPORT);

@ISA    = qw(Exporter);
@EXPORT = qw(sum_bam_reads male_female_check );




=head2 sum_bam_reads

  Arg [1]   : bam file name
  Arg [2]   : Chromosome region to count read for
  Arg [3]   : Boolean true/false for verbose output
  Function  : Count the number of reads for each chromosome in a bam file
  Returntype: reference to hash containings count keys by CHR designation
  Exceptions: throw if file does not exist
  Example   : sum_bam_reads( $bam , 16) or sum_bam_reads( $bam , 16, 1);
=cut
sub sum_bam_reads {
#Note: Probably not quickest way to do this. Effectively cat'ing potentially, very
#      large files 22 times might be slow

  my ($file,  $chr,  $verbose) = @_;

  if (! $chr){
    $chr = [ 1..20, 'X', 'Y', "MT"] ;
  }
  else{
    $chr = [$chr];
  }


  my $samtools = "/nfs/1000g-work/G1K/work/bin/samtools/samtools view ";
  my $wcl      = " | wc -l";
  my $ctr      = 0;
  my %totals;
 
   if (! (-e $file) ){
     warn ("$file doesn not exist. Skipping");
     exit;
   }

   foreach my $chk (@$chr){
     print "running: $samtools $file $chk  $wcl\n" if $verbose;
     $ctr = `$samtools $file $chk  $wcl `;
     chomp $ctr;

     $totals{$chk} = $ctr;

     print "$file chr $chk count = $ctr\n" if $verbose;
     $ctr = 0;
   }
  
  foreach my $key (keys %totals){
    print  $key , "   " ,$totals{$key},"\n" if $verbose;
  }
  return \%totals;
}




#Much faster version of above function. Uses Samtools "pileup".
#Only goes through bam file once.
#example usage: my ( $counts, $sex) = male_female_check ($file, $verbose);
#returns reference to hash containing read counts and also male/female
#designation
sub male_female_check {
my $file    = shift;
my $verbose = shift if (@_);
my $sex;
my $ctr = 0;

throw "$file:File does not exist" if (!-e $file);

if (! ($ENV{'PATH'} =~ /samtools/i) ) {
  throw "samtools location not in \$PATH. Must append";
}

my @chr = ( 1..20, 'X', 'Y', "MT") ;
 

my %chrom_reads;
foreach (@chr){
  $chrom_reads{$_} = 0;
}
 
print "Checking designation of sample\n"          if $verbose;

open(BAMCAT,"samtools pileup $file |"); 
  my $i=0; 
   while (<BAMCAT>) { 
      chomp; 
       my @data      = split(' ', $_);
       $chrom_reads { $data[0]}++; 
   } 
   close(BAMCAT); 

if ($verbose){
  foreach my $key ( keys %chrom_reads){
    print  $key,"\t" ,$chrom_reads{$key},"\n"  if( $chrom_reads {$key} > 0) ;
  }
}

$sex = "male"    if ($chrom_reads{Y} > 0);
$sex = "female"  if ($chrom_reads{Y} == 0);


#print " in sub sex = $sex  $chrom_reads{Y}\n";  

return  \%chrom_reads, $sex  ;

}##




