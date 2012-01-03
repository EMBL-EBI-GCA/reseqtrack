#!/sw/arch/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::BamUtils;
use ReseqTrack::Tools::Intersection;

$| = 1;

my $alignment_index;
my $seq_index;
my $root = '/nfs/1000g-archive/vol1/ftp/';
my $count = 0;
my %collection_bas_files;
my %collection_RG;
my %collection_files;

&GetOptions(
  'alignment_index:s' => \$alignment_index,
  'ftp_root:s' => \$root,
  'seq_index=s' => \$seq_index,
    );

=head Example lines of an alignment index
BAM FILE        BAM MD5 BAI FILE        BAI MD5 BAS FILE        BAS MD5
data/HG00096/alignment/HG00096.chrom1.ILLUMINA.bwa.GBR.low_coverage.20100901.bam        47534db41bd2ae82dc0ab2466d88690c        data/HG00096/alignment/HG00096.chrom1.ILLUMINA.bwa.GBR.low_coverage.20100901.bam.bai64186436d65503f9171ddbc430c818dc data/HG00096/alignment/HG00096.chrom1.ILLUMINA.bwa.GBR.low_coverage.20100901.bam.bas    292c0b96b04912473e25cabe717ca6e0
data/HG00096/alignment/HG00096.chrom10.ILLUMINA.bwa.GBR.low_coverage.20100901.bam       8f0027981ad502fd6d8f7d8428b625c2        data/HG00096/alignment/HG00096.chrom10.ILLUMINA.bwa.GBR.low_coverage.20100901.bam.ba65101d088998d3a08fed419b6706f3a7 data/HG00096/alignment/HG00096.chrom10.ILLUMINA.bwa.GBR.low_coverage.20100901.bam.bas   326793a59bfef53772a6e80998755a33
=cut

my $lines = get_lines_from_file($alignment_index);

foreach my $line(@$lines){
  if($count == 0){
    unless($line =~ /BAM FILE/){
      print STDERR "First line appears not to be a header\n";
      print $line."\n";
    }
    $count++;
    next;
  }
  $count++;
  my @values = split /\t/, $line;
  unless(@values == 6){
    print STDERR $count." ".$line."\n";
    print STDERR "has the wrong number of columns\n";
  }
  my $bam_path = $root.$values[0];
  my $bam_md5 = $values[1];
  my $bai_path = $root.$values[2];
  my $bai_md5 = $values[3];
  my $bas_path = $root.$values[4];
  my $bas_md5 = $values[5];
  
  #next if ($bam_path =~ /20101123/ && $bam_path =~ /chrom20/);  # for the BAMs that have been edited by Shane, he did not replace the chrom20 files.
  																# for 20101123 release, only check whole genome BAMs (mapped and unmapped, not chrom20 ones) 
  
  next if ( check_values($bam_path, $bam_md5) == 1 ) ;
  next if ( check_values($bai_path, $bai_md5) == 1 ) ;
  next if ( check_values($bas_path, $bas_md5) == 1 );
  
  my ($collection_name, $sample, $platform, $algorithm, $project, $analysis, $chrom, $date) = get_collection_name_from_file_name($bam_path);
  $collection_files{$collection_name}{$bam_path} = 1;
  $collection_bas_files{$collection_name}{$bas_path} = 1;
}

## Populate hash %collection_RG ##
my %collection_read_cnt_bas;
my %collection_mapped_read_cnt_bas;
my %collection_base_cnt_bas;
foreach my $col (keys %collection_bas_files) {
	
	my @all_bas_files = keys %{$collection_bas_files{$col}};
	
	my $file_cnt_per_collection = keys %{$collection_files{$col}};
	
	if ( $file_cnt_per_collection != 3 && $col =~ /low_coverage/ ) {
		print STDERR "EACH low_coverage collection should have 3 BAM files, $col has $file_cnt_per_collection bam files\n";
		next;
	}
	elsif ( $file_cnt_per_collection != 1 && $col =~ /exome/ ) {
		print STDERR "EACH exome collection should have 1 BAM file, $col has $file_cnt_per_collection bam files\n";
		next;
	}
	
	$collection_read_cnt_bas{$col} = 0;
	$collection_mapped_read_cnt_bas{$col} = 0;
	$collection_base_cnt_bas{$col} = 0;			
				
	foreach my $bas ( @all_bas_files ) {
		
		next if ( $bas =~ /chrom20/ );
		open (BAS, "<", $bas) || throw("Cannot open bas file $bas\n");
		
		while ( <BAS> ) {
			chomp;
			next if ($_ =~ /bam_filename/); #skip the bas header
			my @tmp = split(/\t/, $_);
			my $RG = $tmp[6];
			my $total_bases = $tmp[7];
			my $read_cnt = $tmp[9];
			my $mapped_read_cnt = $tmp[10];
			$collection_RG{$col}{$RG} = 1;
			$collection_read_cnt_bas{$col} += $read_cnt;
			$collection_base_cnt_bas{$col} += $total_bases;
			$collection_mapped_read_cnt_bas{$col} += $mapped_read_cnt;	
		}
	}
}	

=header -- BAS example lines
bam_filename    md5     study   sample  platform        library readgroup       #_total_bases   #_mapped_bases  #_total_reads   #_mapped_reads  #_mapped_reads_paired_in_sequencing     #_mapped_reads_properly_paired    %_of_mismatched_bases   average_quality_of_mapped_bases mean_insert_size        insert_size_sd  median_insert_size      insert_size_median_absolute_deviation   #_duplicate_reads
NA12275.mapped.ILLUMINA.BWA.CEU.exome.20110521  bf719900c578f5d109b9104ed55e49b8        exome   NA12275 ILLUMINA        2864119779      SRR078850       4893178000      4704651533      48931780        48461017  48461017        47794955        0.24    37.01   8222    883076.31       291     75      3840141
NA12275.mapped.ILLUMINA.BWA.CEU.exome.20110521  bf719900c578f5d109b9104ed55e49b8        exome   NA12275 ILLUMINA        2864119779      SRR078851       4939068000      4746269803      49390680        48910112  48910112        48232632        0.24    37.14   8077    874265.61       291     75      3786178
=cut		  

my ($col_rg_seq_in_ref, $col_read_cnt_seq_in_ref, $col_base_cnt_seq_in_ref) = PARSE_INDEX($seq_index);
my %col_rg_seq_in = %$col_rg_seq_in_ref;
my %col_read_cnt_seq_in = %$col_read_cnt_seq_in_ref;
my %col_base_cnt_seq_in = %$col_base_cnt_seq_in_ref;

foreach my $collection ( keys %collection_RG) {
	
	my $rg_count = keys %{$collection_RG{$collection}};

	my ($sample, $plat, $mapper, $analysis) = split(/\./, $collection);
	$plat = "ILLUMINA" if ($plat =~ /illumina/);
	my $rg_cnt_seq_in = keys %{$col_rg_seq_in{$sample}{$plat}{$analysis}};
	my $read_cnt_seq_in =  $col_read_cnt_seq_in{$sample}{$plat}{$analysis};
	my $base_cnt_seq_in = $col_base_cnt_seq_in{$sample}{$plat}{$analysis};
	#print "Processing collection $collection....\n";
	if ( $rg_cnt_seq_in != $rg_count ) {
		
		print "RG cnt inconsistent: $sample collection $collection has $rg_count RG in BAS "; 
		print "$rg_cnt_seq_in RG in sequence index\n";
		
		my @list_based_on_bas = keys %{$collection_RG{$collection}};
		my @list_based_on_seq_index = keys %{$col_rg_seq_in{$sample}{$plat}{$analysis}};
		my $list_set1 = ReseqTrack::Tools::Intersection
                     ->new(
                           -LIST => \@list_based_on_bas,
                          );
        my $list_set2 = ReseqTrack::Tools::Intersection
                     ->new(
                           -LIST => \@list_based_on_seq_index,
                          );                  
        my $unique_to_bas = $list_set1->not($list_set2);
        my $unique_to_index = $list_set2->not($list_set1);
        
		print "Unique to BAS: ";
		print join("\t", @{$unique_to_bas->list}) . "\n";
		
		print "Unique to SEQ: ";
		print join("\t", @{$unique_to_index->list} ) . "\n";
	}	
	## If the RG count is the same, check read counts	
	else {
		if ($collection_read_cnt_bas{$collection} != $read_cnt_seq_in ) {
			if ( $collection =~ /low_coverage/i ) { ## low_coverage BAMs should have unmapped BAMs
				print "Read cnt inconsistent: $sample collection $collection has $collection_read_cnt_bas{$collection} reads in BAS "; 
				print "$read_cnt_seq_in reads in sequence index\n";
			}
			else {	## EXOME BAMs usually don't have unmapped BAMs	
				my $diff_ratio = ( $read_cnt_seq_in - $collection_read_cnt_bas{$collection})/$read_cnt_seq_in;
				if ($diff_ratio > 0.25) {
					print "Missing more than 25% reads: $sample collection $collection has $collection_read_cnt_bas{$collection} reads in BAS "; 
					print "$read_cnt_seq_in reads in sequence index\n";
				}	
			}	
		}
		## if read count is the same, check for total base count
		else {
			if ($collection_base_cnt_bas{$collection} != $base_cnt_seq_in ) {		
				if ( $collection =~ /low_coverage/i ) {
					print "Base cnt inconsistent: $sample collection $collection has $collection_base_cnt_bas{$collection} bases in BAS "; 
					print "$base_cnt_seq_in total bases in sequence index\n";
				}	
			}
			else {
				print "Pass: $sample collection $collection\n";
			}
		}  
	}			
	
}	

### SUB ###
sub check_values{
  my ($path, $md5) = @_;
  my $flag = 0;
  unless(-e $path){
    print STDERR $path." does not exist\n";
  	$flag = 1;
  }
  unless(length($md5) == 32){
    print STDERR "There is a problem with the md5 value for ".$path."\n";
  	$flag = 1;
  }
  return $flag;
}

sub PARSE_INDEX {
=head
		1.  FASTQ_FILE, path to fastq file on ftp site  
        2.  MD5, md5sum of file
        3.  RUN_ID, SRA/ERA run accession       
        4.  STUDY_ID, SRA/ERA study accession   
        5.  STUDY_NAME, Name of stury   
        6.  CENTER_NAME, Submission centre name 
        7.  SUBMISSION_ID, SRA/ERA submission accession 
        8.  SUBMISSION_DATE, Date sequence submitted, YYYY-MM-DAY       
        9.  SAMPLE_ID, SRA/ERA sample accession 
        10. SAMPLE_NAME, Sample name    
        11. POPULATION, Sample population, this is a 3 letter code and it is defined in README.populations     
        12. EXPERIMENT_ID, Experiment accession 
        13. INSTRUMENT_PLATFORM, Type of sequencing machine     
        14. INSTRUMENT_MODEL, Model of sequencing machine       
        15. LIBRARY_NAME, Library name  
        16. RUN_NAME, Name of machine run       
        17. RUN_BLOCK_NAME, Name of machine run sector  
        18. INSERT_SIZE, Submitter specifed insert size 
        19. LIBRARY_LAYOUT, Library layout, this can be either PAIRED or SINGLE 
        20. PAIRED_FASTQ, Name of mate pair file if exists (Runs with failed mates will have 
            a library layout of PAIRED but no paired fastq file)
        21. WITHDRAWN, 0/1 to indicate if the file has been withdrawn, only present if a file has been withdrawn
        22. WITHDRAWN_DATE, date of withdrawal, this should only be defined if a file is 
            withdrawn
        23. COMMENT, comment about reason for withdrawal
        24. READ_COUNT, read  count for the file
        25. BASE_COUNT, basepair count for the file
        26. ANALYSIS_GROUP, the analysis group of the sequence, this reflects sequencing
            strategy. Currently this includes low coverage, high coverage, exon targetted and exome
            to reflect the 2 non low coverage pilot sequencing stratergies and the 2 main project sequencing stratergies used by the 
            1000 genomes project.
=cut
	my ($index) = @_;
	
	my (	%col_RG, 
			%col_read_cnt,
			%col_base_cnt,
		);
			 
	open (INDEX, "<", $index) || throw("Cannot open sequence index file $index\n");
	
	while (<INDEX>) {
		chomp;
		next if ($_ =~ /FASTQ_FILE/);
		my @data = split (/\t/, $_);
		next if ($data[20] == 1); # jump over fastq files that have been withdrawn
		
		$data[25] = "low_coverage" if ($data[25] =~ /low coverage/);
		$data[25] = "high_coverage" if ($data[25] =~ /high coverage/);
		$data[25] = "exon_targetted" if ($data[25] =~ /exon targetted/);
		$data[12] = "SOLID" if ( $data[12] =~ /ABI_SOLID/ );
		
		# first key is sample name, second key is platform, third key is analysis group, fourth key is RG id
		$col_RG{$data[9]}{$data[12]}{$data[25]}{$data[2]} = 1;		
		
		#### first key is sample name, second key is platform, third key is analysis group, value is the total number of read cnt for the collection
		$col_read_cnt{$data[9]}{$data[12]}{$data[25]} = 0 if (! $col_read_cnt{$data[9]}{$data[12]}{$data[25]});
		if ($data[23] =~ /^\d+$/) {
			$col_read_cnt{$data[9]}{$data[12]}{$data[25]} = $col_read_cnt{$data[9]}{$data[12]}{$data[25]} + $data[23]; 
		}
		
		#### first key is sample name, second key is platform, third key is analysis group, value is the total base count for the collection
		$col_base_cnt{$data[9]}{$data[12]}{$data[25]} = 0 if (! $col_base_cnt{$data[9]}{$data[12]}{$data[25]});
		if ($data[24] =~ /^\d+$/) {
			$col_base_cnt{$data[9]}{$data[12]}{$data[25]} = $col_base_cnt{$data[9]}{$data[12]}{$data[25]} + $data[24]; 
		}	
	}
	return (\%col_RG, \%col_read_cnt, \%col_base_cnt);
}			

=head
 perl /nfs/1000g-work/G1K/work/zheng/reseqtrack/scripts/process/check_alignment_index_sanity.pl -alignment_index 20101123.alignment.index -ftp_root /nfs/1000g-archive/vol1/ftp/ -seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20101123.sequence.index > 20101123.alignment.index.sanity &

 perl /nfs/1000g-work/G1K/work/zheng/reseqtrack/scripts/process/check_alignment_index_sanity.pl -alignment_index 20101123.alignment.index -ftp_root /nfs/1000g-archive/vol1/ftp/ -seq_index /nfs/1000g-archive/vol1/ftp/sequence.index > 20101123.alignment.index.sanity &
 