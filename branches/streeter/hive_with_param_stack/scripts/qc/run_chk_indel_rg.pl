#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::Tools::Exception;
use File::Basename;
use File::Path;
use Getopt::Long;
use Time::Local;

my $bam;
my $outdir;

&GetOptions(
  'bam=s'     		=> \$bam,
  'outdir=s' 		=> \$outdir,	
);

if ($bam !~ /vol1/ ) {
	throw("Cannot run $bam; only check BAMs that have been placed on the ftp site");
}
	
unless (-e $bam) {
	throw("Bam $bam does not exist");
}

if ($bam =~ /chrom|unmapped/i) {
	goto END;
}

#if ($bam !~ /20111114/ ) {
#	goto END;
#}	

	
$outdir =~ s/\/$//;
mkpath $outdir unless (-e $outdir);
my $outfile = $outdir . "/" . basename($bam)  . ".out";

`/nfs/1000g-work/G1K/work/bin/samtools_dev/samtools/chk_indel_rg $bam > $outfile`;

my $exit = $?>>8;
throw("mv failed\n") if ($exit >=1);

END:

=pod
perl $ZHENG_RT/scripts/qc/run_chk_indel_rg.pl -bam /nfs/1000g-archive/vol1/ftp/data/NA20757/alignment/NA20757.chrom20.ILLUMINA.bwa.TSI.low_coverage.20111114.bam -outdir /nfs/1000g-work/G1K/scratch/bam_release_20111114/chk_indel_rg_output