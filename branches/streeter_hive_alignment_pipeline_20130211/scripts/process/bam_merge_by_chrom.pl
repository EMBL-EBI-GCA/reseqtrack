#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use Getopt::Long;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use Time::Local;
use VertRes::Utils::Sam;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::BamUtils;
use File::Path;
use File::Basename;

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $bam_list,
    $collection_type,
    $collection_name,
    $store,
	%bams_by_col_chrom, #hash of arrays, first key is collection name; second key is chrom, value is a list of split BAMs need to be merged
	$only, #merge bams belong to only a certain chr
	$mirror, #when split bams are on a mirror site
	$verbose,
);
    
my $output_dir = "./";
my $start_time = time();
my $host_name = "1000genomes.ebi.ac.uk";
my $remote = 0;

&GetOptions( 
  'dbhost=s'      		=> \$dbhost,
  'dbname=s'      		=> \$dbname,
  'dbuser=s'      		=> \$dbuser,
  'dbpass=s'      		=> \$dbpass,
  'dbport=s'      		=> \$dbport,
  'bam_list:s'			=> \$bam_list,
  'output_dir:s'		=> \$output_dir,
  'collection_type:s'	=> \$collection_type,
  'collection_name:s'	=> \$collection_name,
  'store!'				=> \$store,
  'only:s'				=> \$only,
  'mirror:s'			=> \$mirror,
  'verbose!'			=> \$verbose,
); 

if ($bam_list && $collection_type) {
	throw("Please provide either a list of bams to merge or a collection_type, not both");
}

my $sam_obj = VertRes::Utils::Sam->new();

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

my $ca = $db->get_CollectionAdaptor;
my $fa = $db->get_FileAdaptor;
my $hist_a = $db->get_HistoryAdaptor;
my $ha = $db->get_HostAdaptor;

my $host = $ha->fetch_by_name($host_name);
if(!$host){
  $host = ReseqTrack::Host->new
      (
       -name => $host_name,
       -remote => $remote
      );
}

$db->dbc->disconnect_when_inactive(1);
	
if($bam_list){
  my $lines = get_lines_from_file($bam_list);
  foreach my $line(@$lines){
		throw($line." does not exist on the current filesystem") unless(-e $line);
 		next if ($line !~ /\.bam$/);
 		if ($line =~ /^chrom/) {          # split bams made by sanger code
			my @bits = split (/\./, $line);
			my $chrom = $bits[0];
			push @{$bams_by_col_chrom{"mock_col"}{$chrom}}, $line;
 		}
 		elsif ( $line =~ /chrom/ ) {		# chrom bam uploaded by a center
 		    my ($sample, $platform, $algorithm, $project, $analysis, $chrom, $date, $pop) = CHECK_AND_PARSE_FILE_NAME($line);
			push @{$bams_by_col_chrom{"mock_col"}{$chrom}}, $line;
 		}
 		else {
 		    warning("input bam has to contain chrom assignment, $line doesn't, cannot process\n");
 		} 		    
  }
  $collection_type = "STANDALONE_TD_BAM";
}
elsif ($collection_type && $collection_name)  {
	throw("Collection type has to be  *_TO_TRANSPOSE") if ($collection_type !~ /TO_TRANSPOSE/);
	my $bams_by_col_chrom_ref = process_col($collection_type, $collection_name);
	%bams_by_col_chrom = %$bams_by_col_chrom_ref;
}		
else {
	throw("Please provide either a list of BAMs to merge or a TO_TRANSPOSE collection_name and type");
}	

foreach my $collection ( keys %bams_by_col_chrom ) {
	my %bams_by_chr = %{$bams_by_col_chrom{$collection}};
	run_merge(\%bams_by_chr, $collection, $output_dir, $collection_type); 
}			

my $time_elapsed = time() - $start_time;

print "Job took $time_elapsed seconds\n";	

##### SUBS #####
sub process_col {
	my ($collection_type, $col_name) = @_;
	my $col = $ca->fetch_by_name_and_type($col_name, $collection_type);
	my %bams_by_chrom1;

	my $others = $col->others;
	throw("No bam files found associating with collection " . $col->name) if (!$others || @$others == 0);
	foreach my $other ( @$others ) {
		my $bam_path = $other->name;
		my $bam_basename = basename($bam_path);
		throw("BAM file before splitting does not exist") unless ( -e $bam_path );
		if ( $mirror ) {  	## this is to use for chrom20 BAMs my process did not split. 
							##They already exist in the db, with ftp path.  However they are also on a mirror site
			$mirror =~ s/\/$//;
			throw("Please specify 'only' when use mirror option") unless ($only);
			my $chr = "chrom" . $only;
			my $split_bam_basename = $bam_basename;
			$split_bam_basename =~ s/mapped/$chr/;
			
			my ($sample, $platform, $algorithm, $project, $analysis, $chrom, $date, $pop) = CHECK_AND_PARSE_FILE_NAME($bam_path);
			
			my $split_bam_on_mirror = $mirror . "/data/$sample/alignment/" . $split_bam_basename;
			throw ("No split BAM $split_bam_on_mirror exists on the mirror site") unless (-e $split_bam_on_mirror);
			
			my $split_bam_objs = $fa->fetch_by_filename($split_bam_basename);	
			print "bams to merge for " . $col->name . " is $split_bam_on_mirror\n" if ($verbose);
			push @{$bams_by_chrom1{$col->name}{$chr}}, $split_bam_on_mirror;
			
		}
		else {	
			for (my $i = 1; $i < 23; $i++) {
			### FIXME, add X and Y!!
				my $chr = "chrom" . $i;
				if ($only) {
					$only = "chrom" . $only;
					next if ($chr ne $only);
				}
				my $split_bam_basename = $chr . "." . $bam_basename;
				my $split_bam_objs = $fa->fetch_by_filename($split_bam_basename);
				if ( !$split_bam_objs  || @$split_bam_objs == 0) {
					throw("No split $chr BAM object exists for BAM $bam_path") if ( !$split_bam_objs  || @$split_bam_objs == 0);
					## FIXME: carry on to split the BAM
				}
				elsif (@$split_bam_objs > 1) {
					throw("More than one $chr split bams exist for bam $bam_path"); 
				}	
				else {
					print "bams to merge for " . $col->name . " is " . $split_bam_objs->[0]->name . "\n" if ($verbose);
					push @{$bams_by_chrom1{$col->name}{$chr}}, $split_bam_objs->[0]->name;
				}
			}
		}
	}	
	return \%bams_by_chrom1;				
}	
	
sub run_merge {
	my ($bams_by_chr_ref, $coll, $out_dir, $col_type) = @_;
	my %bams_by_chr = %$bams_by_chr_ref;
	$out_dir =~ s/\/$//;
	
	mkpath($out_dir) unless (-e $out_dir);
	
	my $type = $col_type;
	$type =~ s/_TO_TRANSPOSE//;
	$type = "TD_" . $type;  ## TD stands for TransformeD
 
	foreach my $chr ( keys %bams_by_chr ) {
		my $merged_bam = $out_dir . "/" . $coll . "_" . $chr . ".bam";
		if ( $sam_obj->merge($merged_bam, @{$bams_by_chr{$chr}}) ) {
			print "Merged bam is $merged_bam\n";
			store_file($merged_bam, $type) if ($store);
			`samtools index $merged_bam`; #create index file, needed for variant calling
			my $exit = $?>>8;
			throw("Cannot index bam file $merged_bam\n") if ($exit >=1);
		}
		else {
			print "Run merge failed\n";
		}
	}
	return 1;
}	

sub store_file {
	my ($file_path, $file_type) = @_;
	my $basename = basename($file_path);
	my $existing_fo = $fa->fetch_by_filename($basename);
	
	if (!$existing_fo || @$existing_fo == 0) {	
		my $fos = create_objects_from_path_list([$file_path], $file_type, $host); #$files is a reference to an array of file objects (in this case, only one element)
		my $fo = $fos->[0]; #to get the first and only file object
	
		if ($store) {
			$fa->store($fo);	
			my $history = ReseqTrack::History->new(
				-other_id 	=> $fo->dbID,
				-table_name => 'file',
				-comment 	=> "merge bams of the same chromosome to form a transformed bam", 
			);
			$hist_a->store($history);
			$fo->history($history);
		}
	}
	else {
		my $new_fo = ReseqTrack::File->new
 	   			(
	      		  -adaptor => $fa,
			      -dbID => $existing_fo->[0]->dbID,
			      -name => $file_path,
#			      -md5 => "NULL",
			      -host => $host,
			      -type => $file_type,
			   );
			
			my $history_ref = $existing_fo->[0]->history;
			my $history;
			
			if (!$history_ref || @$history_ref == 0) {     	  
				$history = ReseqTrack::History->new(
				-other_id => $existing_fo->[0]->dbID,
				-table_name => 'file',
				-comment => "Update transformed BAM file", 
				);
				$new_fo->history($history);
			}
			else {
				my $comment = calculate_comment($existing_fo->[0], $new_fo); 
				
				if (!$comment) {
					$comment = "fiddling tranformed bam files, no comment\n";
				}	
				$history = ReseqTrack::History->new(
					-other_id => $existing_fo->[0]->dbID,
					-table_name => 'file',
					-comment => $comment,
				);
				$new_fo->history($history);	
			}
			
			if ($store) {
				$fa->update($new_fo, 1, 1); # the second 1 is allow change name		
			}
	}
}	

#####################################################################################################################################################
=pod

=head1 SYNOPSIS
If -collection_name not given, all collections with the tyoe defined by -collection_type would be processed

perl $ZHENG_RP/bin/bam_merge_by_chrom.pl $WRITE_DB_ARGS -dbname zheng_var_call -collection_type PHASE1_BAM_TO_TRANSPOSE -collection_name low_coverage_PUR_SOLID -output_dir /nfs/1000g-work/G1K/scratch/zheng/transpose_bam/test_merge/ -store & 

perl /nfs/1000g-work/G1K/work/zheng/reseq-personal/zheng/bin/bam_merge_by_chrom.pl -dbhost mysql-g1kdcc-public -dbname zheng_var_call -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -output_dir /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20 -collection_type PHASE1_BAM_TO_TRANSPOSE  -mirror /nfs/1000g-work/G1K/work/phase1_alignments -only 20 -store -collection_name low_coverage_PUR_SOLID