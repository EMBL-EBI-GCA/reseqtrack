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
    $input_bam,
    $pretend,
    $only,
    $run,
    $mirror,  # use this when unsplit BAMs are copied to a mirror site on G1K, not on ftp site, something like /nfs/1000g-work/G1K/work/phase1_alignments
);
    
my $output_root = "./";
my $start_time = time();
my $host_name = "1000genomes.ebi.ac.uk";
my $remote = 0;
my $keep_file_type = 0; 		#if want to keep the file type unchanged
my $phase2 = 0; 				# this is to indicate that the split is for phase 2 Baylor EXOME BAMs
my $edit_split_bam_name = 0; 	#this is to indicate whether you want to use conventional bam name for split bams.  
								#unsplit BAM need to fit naming convention to begin with 
my $make_bai = 0; 				#to indicate if bai file is needed for the split bam

&GetOptions( 
  'dbhost=s'      		=> \$dbhost,
  'dbname=s'      		=> \$dbname,
  'dbuser=s'      		=> \$dbuser,
  'dbpass=s'      		=> \$dbpass,
  'dbport=s'      		=> \$dbport,
  'bam_list=s'			=> \$bam_list,
  'bam=s'				=> \$input_bam,
  'output_root=s'		=> \$output_root,
  'pretend!'			=> \$pretend,
  'only=s'				=> \$only,
  'run!'				=> \$run,
  'mirror=s'			=> \$mirror,
  'remote=i'			=> \$remote,
  'keep_file_type!'		=> \$keep_file_type,
  'edit_bam_name!'		=> \$edit_split_bam_name,
  'is_phase2!'			=> \$phase2,
  'host:s'				=> \$host_name,
  'make_bai!'			=> \$make_bai,
); 

$output_root =~ s/\/$//;

my $sam_obj = VertRes::Utils::Sam->new();

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

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

my @input_bams;
if($bam_list){
  my $lines = get_lines_from_file($bam_list);
  foreach my $line(@$lines){
		throw($line." does not exist on the current filesystem") unless(-e $line);
		throw("This is already a chrom BAM, skip $input_bam") if ( $input_bam =~ /chrom/ );
 		next if ($line !~ /\.bam$/);
	    push(@input_bams, $line);
  }
}
elsif ($input_bam) {
	throw("BAM $input_bam doesn't exist!") unless (-e $input_bam);
	if ( $input_bam =~ /chrom/ || $input_bam =~ /unmapped/ ) {
		goto END;
		warning("BAM not split by chrom; this is a chrom BAM or a unmapped BAM, skip $input_bam");
	}	
	
	if ( $phase2 ) {
		if ( $input_bam !~ /SOLID\.bfast/i ) { # to skip all non-baylor exome bams; they will be labeled as complete in EC table
			goto END;
			warning("BAM not split. It is not Baylor exome bams. BAM is $input_bam");
		}	
		else { # For baylor BAM split
			if ($input_bam !~ /vol1/ || $input_bam !~ /exome\.20111114/ ) {
				throw("BAM not split. To be split, phase 2 BAM needs to be passed DCC QA and uploaded to the ftp site. BAM is $input_bam");
			}
		}	
	}
	push @input_bams, $input_bam;
}
else {
	throw("Please provide either a path to a BAM file or a file containing a list of bam files");
}		

foreach my $bam ( @input_bams ) {
	my $file_obj = $fa->fetch_by_name($bam);
	my $split_file_type;
	if ($keep_file_type ) {
		$split_file_type = $file_obj->type;
	}
	else {	
		$split_file_type = "SPLIT_" . $file_obj->type;
	}
	my ($sample, $platform, $algorithm, $project, $analysis, $chrom, $date, $pop) = CHECK_AND_PARSE_FILE_NAME($bam);
	my $grp = $analysis . "_" . "$pop" . "_" . $platform;
	my $output_dir = $output_root . "/$grp/";
	mkpath($output_dir) unless (-e $output_dir);
	my $config_hash = {
		'output_dir' => $output_dir,
		'pretend'	 => $pretend,
		'only'		=> $only,
	};	##hash_ref
	
	my @bams;
	if ($mirror) {
		$mirror =~ s/\/$//;
		my $bam_basename = basename($bam);
		my $bam_on_mirror = $mirror . "/data/$sample/alignment/" . $bam_basename;
		### FIXME, use exome_alignment if necessary
		throw("BAM $bam_on_mirror does not exist on the mirrored site") unless (-e $bam_on_mirror); 
		@bams = $sam_obj->split_bam_by_sequence($bam_on_mirror, %$config_hash);	
	}
	else {		
		@bams = $sam_obj->split_bam_by_sequence($bam, %$config_hash);
	}
	## without any specifications, it generates umapped, nonchrom, X, Y, MT and all autosomes
	## FIXME, if there is chrom file exist in the designated location, then no split will happen! need to change this behavior??

	foreach my $split_bam ( @bams ) {
		if ($split_bam =~ /nonchrom/ && $phase2) {
			`rm $split_bam` if (-e $split_bam);
			my $exit3 = $?>>8;
			throw("rm failed\n") if ($exit3 >=1);
			next;
		}	
		print "split bam is $split_bam\n";
		
		my $edited_split_bam_path;
		my $bai;
		my $bai_type;
		if ( $edit_split_bam_name ) {
			my $dir = dirname($split_bam);
			my $bname = basename($split_bam);
			my @bits = split(/\./, $bname );
			my $chrom = shift @bits;
			my $edited_split_bam = join(".", @bits);
			$edited_split_bam =~ s/mapped/$chrom/;
			$edited_split_bam_path = $dir . "/$edited_split_bam";
			print "edited split bam is $edited_split_bam_path\n";	
			`mv $split_bam 	$edited_split_bam_path` if ($run);
			my $exit = $?>>8;
			throw("mv failed\n") if ($exit >=1);
			
			if ($make_bai) {
				$bai = make_bai($edited_split_bam_path);
				my $exit2 = $?>>8;
				throw("make_bai failed\n") if ($exit2 >=1);
				$bai_type = $split_file_type;
				$bai_type =~ s/BAM/BAI/;
			}	
		}
		
		if ( $edit_split_bam_name ) {
			store_file($edited_split_bam_path, $split_file_type) unless ($pretend);
			store_file($bai, $bai_type) if (!$pretend && $make_bai);
		}
		else {
			store_file($split_bam, $split_file_type) unless ($pretend);
		}
	}	
}

END: 

my $lasped_time = time() - $start_time;
print "job took $lasped_time seconds\n";

###### SUBS ######	
sub make_bai {
	my ($bam_path) = @_;
	my $output_bai = $bam_path . ".bai";
	`/nfs/1000g-work/G1K/work/bin/samtools/samtools index $bam_path $output_bai`;
	return $output_bai;
}
	 	
sub store_file {
	my ($file_path, $file_type) = @_;
	my $basename = basename($file_path);
	my $existing_fo = $fa->fetch_by_filename($basename);
	
	if (!$existing_fo || @$existing_fo == 0) {	
		my $fos = create_objects_from_path_list([$file_path], $file_type, $host); #$files is a reference to an array of file objects (in this case, only one element)
		my $fo = $fos->[0]; #to get the first and only file object
	
		if ($run) {
			my $md5 = run_md5($file_path);
			$fo->md5($md5);
			$fa->store($fo);	
			my $history = ReseqTrack::History->new(
				-other_id 	=> $fo->dbID,
				-table_name => 'file',
				-comment 	=> "split a whole genome BAM by chromosome", 
			);
			$hist_a->store($history);
			$fo->history($history);
		}
	}
	else {
		my $md5_b = run_md5($file_path);
		my $new_fo = ReseqTrack::File->new
 	   			(
	      		  -adaptor => $fa,
			      -dbID => $existing_fo->[0]->dbID,
			      -name => $file_path,
			      -md5 => $md5_b,
			      -host => $host,
			      -type => $file_type,
			   );
			
			my $history_ref = $existing_fo->[0]->history;
			my $history;
			
			if (!$history_ref || @$history_ref == 0) {     	  
				$history = ReseqTrack::History->new(
				-other_id => $existing_fo->[0]->dbID,
				-table_name => 'file',
				-comment => "Update split BAM file", 
				);
				$new_fo->history($history);
			}
			else {
				my $comment = calculate_comment($existing_fo->[0], $new_fo); 
				
				if (!$comment) {
					$comment = "fiddling split bam files, no comment\n";
				}	
				$history = ReseqTrack::History->new(
					-other_id => $existing_fo->[0]->dbID,
					-table_name => 'file',
					-comment => $comment,
				);
				$new_fo->history($history);	
			}
			
			if ($run) {
				$fa->update($new_fo, 1, 1); # the second 1 is allow change name		
			}
	}
}	

####################################################################################################################################################
=pod

=head1 NAME
bam_split_by_chrom.pl

=head1 SYNOPSIS
This script takes a bam path or a list of bam paths, then split them by chromosome. Output split bams are stored in a folder named by collection 
(analysis_pop_platform); the files are stored in database as type SPLIT_*

=head1 EXAMPLE
perl /nfs/1000g-work/G1K/work/zheng/reseq-personal/zheng/bin/bam_split_by_chrom.pl $WRITE_DB_ARGS -dbname zheng_var_call -bam /nfs/1000g-archive/vol1/ftp/phase1/data/HG01105/alignment/HG01105.mapped.SOLID.bfast.PUR.low_coverage.20101123.bam -output_root /nfs/1000g-work/G1K/scratch/zheng/transpose_bam/PUR_SOLID -mirror /nfs/1000g-work/G1K/work/phase1_alignments -pretend

To split Baylor phase2 exome bams to get chrom20 and chrom11

perl /nfs/1000g-work/G1K/work/zheng/reseq-personal/zheng/bin/bam_split_by_chrom.pl $WRITE_DB_ARGS -dbname zheng_g1k_tracking -bam /nfs/1000g-archive/vol1/ftp/data/NA20529/exome_alignment/NA20529.mapped.SOLID.bfast.TSI.exome.20111114.bam -output_root /nfs/1000g-work/G1K/drop/g1k-drop-bcm -only 11 -keep_file_type -edit_bam_name -is_phase2 -host baylor 