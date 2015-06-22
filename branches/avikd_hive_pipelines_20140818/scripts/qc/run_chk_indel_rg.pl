#!/usr/bin/env perl -w

use strict;

use ReseqTrack::Tools::Exception qw(throw);
use Getopt::Long;
use ReseqTrack::Tools::QC::ChkIndelsBAM qw(run_program);
use ReseqTrack::DBSQL::DBAdaptor;

my $bam;
my $outdir;
my $program;
my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;
my $help;
my $store;

&GetOptions(
  'dbhost=s'     		=> \$dbhost,
  'dbname=s'     		=> \$dbname,
  'dbuser=s'     		=> \$dbuser,
  'dbpass=s'     		=> \$dbpass,
  'dbport=s'     		=> \$dbport,
  'help!'		 		=> \$help,
  'store!'				=> \$store,
  'program=s'			=> \$program,
  'bam=s'	     		=> \$bam,
  'outdir=s' 			=> \$outdir,	
);

$program = '/nfs/1000g-work/G1K/work/bin/samtools_dev/samtools/chk_indel_rg' if (!$program);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
	-host => $dbhost,
	-user => $dbuser,
	-port => $dbport,
	-dbname => $dbname,
	-pass => $dbpass,
);
			
if ($bam !~ /vol1/ ) {
	throw("Cannot run $bam; only check BAMs that have been placed on the ftp site");
}

if ($bam =~ /\.mapped/i) {
	
	my $chk_indels = ReseqTrack::Tools::QC::ChkIndelsBAM->new (
		-db				=> $db,
		-program 		=> $program,
		-input_files 	=> [$bam],
		-working_dir	=> $outdir,
		-store			=> $store,
	);	  
	
	$chk_indels->run;

	my $outfile = $chk_indels->output_files->[0];
	print "out file is $outfile\n";
}
else {
	print "Skip $bam; only check mapped BAMs\n";
}	


=pod

=Head1 NAME

reseqtrack/scripts/qc/run_chk_indel_rg.pl 

=Head1 Synopsis

Excessive amount of short insertions/deletions in a BAM file is a good indication of poor data quality. 
Heng Li developed a script chi_indel_rg.pl to calculate indel content of BAM files. 
A RunProgram style module ChkIndelsBAM can be used to run the script ReseqTrack::Tools::QC::ChkIndelsBAM

This is a caller script developed for the 1KG project.

=Head1 OPTIONS

	-bam		bam file path, have to be on the ftp site; has to be mapped BAMs
	-outdir
	-program
	-store
	-dbname
	-dbhost  
	-dbuser
	-dbpass
	-dbport

=Head1 EXAMPLE COMMAND LINE

perl $ZHENG_RT/scripts/qc/run_chk_indel_rg.pl \
-bam /nfs/1000g-archive/vol1/ftp/data/HG01377/alignment/HG01377.mapped.ILLUMINA.bwa.CLM.low_coverage.20120522.bam \
-outdir /nfs/1000g-work/G1K/scratch/bam_release_20111114/chk_indel_rg_output
$WRITE_DB_ARGS \
-dbname zheng_g1k_copy \
-store

=Head1 OUTPUT

input_file.out in the outdir specified. The header of the out put file is:
file name, RG-ID, #<=6bp-ins, #<=6bp-del, #>6bp-ins and #>6bp-del, 
All results are also stored in the database statistic table.
For 1KG BAMs, anything with ratio of col3/col4 or col4/col3 > 5 is classified as bad. 
The results can be queried using script ~/scripts/qc/query_chk_indel_results.pl

=cut