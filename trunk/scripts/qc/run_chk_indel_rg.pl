#!/usr/bin/env perl -w

use strict;

use ReseqTrack::Tools::Exception qw(throw);
use Getopt::Long;
use ReseqTrack::Tools::QC::ChkIndelsBAM qw(run_program);

my $bam;
my $outdir;
my $program;

&GetOptions(
	'program=s'			=> \$program,
  	'bam=s'     		=> \$bam,
  	'outdir=s' 			=> \$outdir,	
);

$program = '/nfs/1000g-work/G1K/work/bin/samtools_dev/samtools/chk_indel_rg' if (!$program);

if ($bam !~ /vol1/ ) {
	throw("Cannot run $bam; only check BAMs that have been placed on the ftp site");
}

if ($bam =~ /\.mapped/i) {
	
	my $chk_indels = ReseqTrack::Tools::QC::ChkIndelsBAM->new (
		-program 		=> $program,
		-input_files 	=> [$bam],
		-working_dir	=>	$outdir,
	);	  
	
	$chk_indels->run;

	my $outfile = $chk_indels->output_files->[0];
	print "out file is $outfile\n";
}
else {
	print "Skip $bam; only check mapped BAMs\n";
}	


=pod
perl $ZHENG_RT/scripts/qc/run_chk_indel_rg.v2.pl -bam /nfs/1000g-archive/vol1/ftp/data/NA20757/alignment/NA20757.mapped.ILLUMINA.bwa.TSI.low_coverage.20130422.bam -outdir /nfs/1000g-work/G1K/scratch/bam_release_20111114/chk_indel_rg_output