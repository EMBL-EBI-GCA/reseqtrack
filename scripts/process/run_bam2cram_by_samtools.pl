#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use IPC::System::Simple qw(system);
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::Exception;
use File::Basename qw(basename);
use File::Path qw(mkpath);
use ReseqTrack::Tools::RunSamtools qw(run_bam_to_cram run_index);

my (
	$bam,
	$reference,
	$out_dir,
);

&GetOptions(
        'bam=s'       => \$bam,
        'reference=s'	=> \$reference,
        'out_dir=s'				=> \$out_dir,
);

my @name_bits = split(/\./, basename($bam));
my $alignment_dir_name;
if (basename($bam) =~ /high_cov/) {
	$alignment_dir_name = "high_cov_alignment";
}
elsif (basename($bam) =~ /exome/) {
	$alignment_dir_name = "exome_alignment";
}	
else {
	$alignment_dir_name = "alignment";
}	

my $final_dir = $out_dir . "/" . $name_bits[0] . "/" . $alignment_dir_name . "/";

mkpath($final_dir) unless (-e $final_dir);

#### reference need to be ALT or PRIMARY, depends on the input BAM files

my $bam_to_cram = ReseqTrack::Tools::RunSamtools->new(
	-input_files	=> [$bam],
	-reference		=> $reference,
	-working_dir	=> $final_dir,
	);
	
$bam_to_cram->run_bam_to_cram;

my $cram = $bam_to_cram->output_files->[0];

my $create_crai = ReseqTrack::Tools::RunSamtools->new(
	-input_files	=> [$cram],
	-working_dir	=> $final_dir,
	);
	
$create_crai->run_index;

my $crai = $create_crai->output_files->[0];

print "cram is $cram, index is $crai\n";

=pod

perl /nfs/1000g-work/G1K/work/zheng/reseqtrack_branch_zheng/align_output_cram/scripts/process/run_bam2cram_by_samtools.pl \
-bam /nfs/1000g-work/G1K/scratch/zheng/map_by_alt_bwamem_biobambam_to_h38_test/SRS000625/alignment/NA12347.alt_bwamem_bbb_hs38DH.20130502.low_coverage.20150302.bam \
-reference /nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
-out_dir ~/tmp	