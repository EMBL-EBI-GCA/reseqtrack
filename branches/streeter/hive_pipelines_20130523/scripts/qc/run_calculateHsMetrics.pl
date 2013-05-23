#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::DBSQL::DBAdaptor; 
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::Exception;
use Getopt::Long;
use File::Basename;

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $col,
    $out_dir,
    $flanking_50bp,
    $out,
    $col_type,
    $bait,
);    

&GetOptions(
  'dbhost=s'     	=> \$dbhost,
  'dbname=s'     	=> \$dbname,
  'dbuser=s'     	=> \$dbuser,
  'dbpass=s'     	=> \$dbpass,
  'dbport=s'     	=> \$dbport,
  'collection=s'	=> \$col,
  'type=s'			=> \$col_type,
  'out=s'			=> \$out_dir,
  'flanking!'		=> \$flanking_50bp,
  'bait=s'			=> \$bait,
);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );
    
my $ca = $db->get_CollectionAdaptor;
my $fa = $db->get_FileAdaptor;

my $col_obj = $ca->fetch_by_name_and_type($col, $col_type);
if (!$col_obj) {
     throw("No collection found for $col, type $col_type");
}
      
$db->dbc->disconnect_when_inactive(1);    

$out_dir =~ s/\/$//g;

my $others = $col_obj->others;
if (! $others || @$others == 0) {
	throw("No files are found in collection $col, type $col_type");
}

foreach my $other ( @$others ) {
	my $bam = $other->name;

	my $basename = basename($bam);
	
	next if ($bam =~ /chrom|unmapped/i );
	
	if ( $bam !~ /vol1/ ) {
		throw("File $bam not on ftp site yet");
	}		
	
	unless (-e $bam) {
		throw("File $bam does not exist\n");
	}		

	if ($flanking_50bp ) {
		$out = $out_dir . "/" . $basename . ".flanking_50bp.stats";
	}
	else {
		$out = $out_dir . "/" . $basename . ".stats";
	}					
	
	my $command = "/usr/bin/java -jar /nfs/1000g-work/G1K/work/bin/picard/CalculateHsMetrics.jar ";
	$command .= "TARGET_INTERVALS=$bait ";
	$command .= "BAIT_INTERVALS=$bait ";
	$command .=  "INPUT=$bam OUTPUT=$out VALIDATION_STRINGENCY=SILENT VERBOSITY=ERROR ";

	print "$command\n";
	
=head
	if ($flanking_50bp ) {
		if ($bam !~ /mosaik/i) {
			$out = $out_dir . "/" . $basename . ".flanking_50bp.stats";
			#$command = $command . "BAIT_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/20110426_exome_add50bp.consensus.fake_anno.sam ";
			#$command = $command . "TARGET_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/20110426_exome_add50bp.consensus.fake_anno.sam ";
			$command = $command . "BAIT_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/reference/20110426_exome_add50bp.consensus.fake_anno.phase2.sam ";
			$command = $command . "TARGET_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/reference/20110426_exome_add50bp.consensus.fake_anno.phase2.sam ";
		}
		else {
			$out = $out_dir . "/" . $basename . ".flanking_50bp.stats";
			#$command = $command . "BAIT_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/20110426_exome_add50bp.consensus.fake_anno.mosaik.sam ";
			#$command = $command . "TARGET_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/20110426_exome_add50bp.consensus.fake_anno.mosaik.sam ";
			$command = $command . "BAIT_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/reference/20110426_exome_add50bp.consensus.fake_anno.mosaik.phase2.sam ";
			$command = $command . "TARGET_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/reference/20110426_exome_add50bp.consensus.fake_anno.mosaik.phase2.sam ";	
		}
	}
	else {
		if ($bam !~ /mosaik/i){	
			$out = $out_dir . "/" . $basename . ".stats";
	#		$command = $command . "BAIT_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/20110225.exome.consensus.annotation.sam ";
	#		$command = $command . "TARGET_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/20110225.exome.consensus.annotation.sam ";
			$command = $command . "BAIT_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/reference/20110225.exome.consensus.annotation.phase2.sam ";
			$command = $command . "TARGET_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/reference/20110225.exome.consensus.annotation.phase2.sam ";	
		}
		else {
			$out = $out_dir . "/" . $basename . ".stats";
			#$command = $command . "BAIT_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/20110225.exome.consensus.annotation.mosaik.sam ";
			#$command = $command . "TARGET_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/20110225.exome.consensus.annotation.mosaik.sam ";
			$command = $command . "BAIT_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/reference/20110225.exome.consensus.annotation.mosaik.phase2.sam ";
			$command = $command . "TARGET_INTERVALS=/nfs/1000g-work/G1K/work/zheng/exome_stats/reference/20110225.exome.consensus.annotation.mosaik.phase2.sam ";
		}	
	}	
=cut
	
	`$command`;
	
	my $exit = $?>>8;
	throw("calculateHsMetrics failed\n") if ($exit >=1);

}	

=pod
perl /nfs/1000g-work/G1K/work/zheng/reseqtrack/scripts/qc/run_calculateHsMetrics.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-out /nfs/1000g-work/G1K/work/zheng/bam_release.20111114/calculateHsMetrics/results \
-bait /nfs/1000g-work/G1K/work/zheng/exome_stats/20110225.exome.consensus.annotation.sam \
-collection NA20801.ILLUMINA.bwa.exome \
-type EXOME_BAM &