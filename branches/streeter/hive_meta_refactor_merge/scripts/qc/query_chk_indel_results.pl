#!/usr/bin/env perl -w

use strict;

use ReseqTrack::Tools::Exception qw(throw);
use Getopt::Long;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(get_lines_from_file);
use ReseqTrack::Tools::BamUtils qw(CHECK_AND_PARSE_FILE_NAME);

my $bam;
my $bam_list;
my $seq_index;
my $outdir;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;
my $help;

&GetOptions(
  'dbhost=s'     		=> \$dbhost,
  'dbname=s'     		=> \$dbname,
  'dbuser=s'     		=> \$dbuser,
  'dbpass=s'     		=> \$dbpass,
  'dbport=s'     		=> \$dbport,
  'help!'		 		=> \$help,
  'bam=s'     			=> \$bam,
  'bam_list:s'			=> \$bam_list,
  'seq_index=s'			=> \$seq_index,
  'outdir=s' 			=> \$outdir,	
);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
	-host => $dbhost,
	-user => $dbuser,
	-port => $dbport,
	-dbname => $dbname,
	-pass => $dbpass,
);

my $seq_index_hash = parse_index($seq_index);
=head
	foreach my $k (keys %$seq_index_hash) {
		print "$k\n";
	}	
=cut

my $rmi_a = $db->get_RunMetaInfoAdaptor;
my $stats_a = $db->get_StatisticsAdaptor;

if ($bam_list ) {
	my $bams = get_lines_from_file($bam_list);
	foreach my $a_bam ( @$bams ) {
		process_results($a_bam);
	}	
}
elsif ( $bam ) {
	process_results($bam);
}	

###### SUBS ######
sub process_results {
	my ($the_bam) = @_;	
	my ($sample, $platform, $algorithm, $project, $analysis, $chrom, $date) = CHECK_AND_PARSE_FILE_NAME($the_bam);

	#print "$sample\t$analysis\n";
	 
	my @run_ids = @{$seq_index_hash->{$sample}->{$analysis}};
	
	my $flag = 0;
	foreach my $run_id ( @run_ids ) {
		print "run id $run_id\t";
		my $rmi_obj = $rmi_a->fetch_by_run_id($run_id);	
		#print "rmi id is " . $rmi_obj->dbID . "\n";
		my $stats_objs = $stats_a->fetch_by_other_id_and_table_name($rmi_obj->dbID, 'run_meta_info');
		if ( ! $stats_objs || @$stats_objs == 0 ) {
			print "no check indel results\n";
			$flag = 1;
			next;
		}	
		my $sh_insert;
		my $sh_del;
		foreach my $stats_obj (@$stats_objs) {
			$sh_insert =  $stats_obj->attribute_value if ($stats_obj->attribute_name eq "sh_insert");
			$sh_del  =  $stats_obj->attribute_value if ($stats_obj->attribute_name eq "sh_del");
		}
		my $ratio1 = $sh_insert/$sh_del;
		my $ratio2 = $sh_del/$sh_insert;
		
		print "$ratio1\t$ratio2\n";
		
		if ($ratio1 > 5 || $ratio2 > 5 ) {
			print "Run $run_id has problem!\n";
			$flag = 1;
		}	
	}		
	
	if ($flag == 0 ) {
		print "BAM file $the_bam passes!\n";
	}
	else {
		print "BAM file $the_bam failed!\n";
	}	
	return 1;
}
					
sub parse_index {
	my ($seq_i) = @_;
	my $lines = get_lines_from_file($seq_i);
	my %sample_analysisGrp_to_runID;
	foreach my $line ( @$lines ) {
		next if $line =~ /FASTQ_FILE/;
		my @data = split(/\t/, $line);
		my $sample_name = $data[9];
		my $run_id = $data[2];
		my $analysis_grp = $data[25];
		$analysis_grp =~ s/ /_/;
		next if ($data[20]==1);	
		push @{$sample_analysisGrp_to_runID{$sample_name}{$analysis_grp}}, $run_id;	
	}
	return \%sample_analysisGrp_to_runID;
}			
			
			
=pod
perl $ZHENG_RT/scripts/qc/query_chk_indel_results.pl $WRITE_DB_ARGS \
-dbname zheng_g1k_copy \
-bam HG01377.mapped.ILLUMINA.bwa.CLM.low_coverage.20120522.bam \
-seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20130502.analysis.sequence.index

OR using a list of BAMs as input:

perl $ZHENG_RT/scripts/qc/query_chk_indel_results.pl $WRITE_DB_ARGS \
-dbname zheng_g1k_copy \
-bam_list test \
-seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20130502.analysis.sequence.index