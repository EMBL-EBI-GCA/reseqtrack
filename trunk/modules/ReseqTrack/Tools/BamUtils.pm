package ReseqTrack::Tools::BamUtils;
use strict;
use warnings;

use Exporter;
use ReseqTrack::File;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use File::Copy;
use File::Basename;
use File::Find ();
use Data::Dumper;

use vars qw (@ISA  @EXPORT);

@ISA    = qw(Exporter);
@EXPORT = qw(sum_bam_reads male_female_check move_bam_to_trash
             CHECK_AND_PARSE_FILE_NAME get_collection_name_from_file_name check_bas
             get_bam_run_id_info);


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

sub  move_bam_to_trash {
	my ($db2, $file, $actual_file_path, $run) = @_;
	my $full_name = $file->name;
	my $filen = basename($full_name);

	##### If a BAM file, remove it from collection
	
	#if ($file->type eq "BAM" || $file->type eq "NCBI_BAM") {
	if ($file->type =~ /BAM/) {
		my ($collection_name) = get_collection_name_from_file_name($full_name);
		
		if ($collection_name) {
			my $ca = $db2->get_CollectionAdaptor;
			my $exist_collection = $ca->fetch_by_name_and_type($collection_name, $file->type);
			my @bams_to_remove;
			push @bams_to_remove, $file;
			$ca->remove_others($exist_collection, \@bams_to_remove) if ($exist_collection && $run);
		}	
	}	

	##### Move the file to reject bin
	my $reject_dir = "/nfs/1000g-work/G1K/drop/reject/";
		
	unless (-e $reject_dir) {
		mkpath($reject_dir);
	}
				
	my $new_file_path = $reject_dir . $filen;
	
	if ( -e  $full_name ) {
		print "old file path is $full_name and file path in reject bin is $new_file_path\n";					
		`mv $full_name $new_file_path` if ($run);	 	
		my $exit = $?>>8;
		throw("mv failed\n") if ($exit >=1);
	}
	else {
		print "old file path is $full_name, physically locates at $actual_file_path, file path in reject bin is $new_file_path\n";
		`mv $actual_file_path $new_file_path` if ($run);
                my $exit = $?>>8;
                throw("mv failed\n") if ($exit >=1);
        }
	
	##### Make the change in the database, change type to REJECTED_BAM, BAS, or BAI, host_id will be 1, history will be changed	
	my $ha = $db2->get_HostAdaptor;
	my $fa = $db2->get_FileAdaptor;
	my $new_host = $ha->fetch_by_name("1000genomes.ebi.ac.uk");
	my $new_type = "REJECTED_" . $file->type;				
	my $new_file_object = ReseqTrack::File->new
 	   (
	      -adaptor => $fa,
	      -dbID => $file->dbID,
	      -name => $new_file_path,
	      -md5 => $file->md5,
	      -host => $new_host,
	      -type => $new_type,
	      -size => $file->size,
	      -created => $file->created
	   );
			     
	my $history_ref = $file->history;
	my $history;
	
	if (!$history_ref || @$history_ref == 0) {     	  
		$history = ReseqTrack::History->new(
			-other_id => $file->dbID,
			-table_name => 'file',
			-comment => "Move file to reject bin while keep md5 unchanged", 
		);
		$new_file_object->history($history);
	}
	else {
		my $comment = ReseqTrack::Tools::FileUtils::calculate_comment($file, $new_file_object); 
		$history = ReseqTrack::History->new(
			-other_id => $file->dbID,
			-table_name => 'file',
			-comment => $comment,
		);
		$new_file_object->history($history);
	}
	
	$fa->update($new_file_object, 1) if($run);
	            	
	return 1;
}

sub get_collection_name_from_file_name {
	my ($file_path) = @_;
	
	my ($sample, $platform, $algorithm, $project, $analysis, $chrom, $date) = CHECK_AND_PARSE_FILE_NAME($file_path);
			
	my $collection_name;
	if ($sample && $platform && $algorithm && $project) {
		$collection_name = join ('.', $sample, $platform, $algorithm, $project); ### date does not have to be the same for files in the same collection
	}
	elsif ($sample && $platform && $algorithm && $analysis) {
		$collection_name = join ('.', $sample, $platform, $algorithm, $analysis); 
	}
	
	return ($collection_name, $sample, $platform, $algorithm, $project, $analysis, $chrom, $date);			
}


#	EXAMPLES of BAM file names in pilots:
#	NA19240.SLX.maq.SRP000031.2009_08.bam
#	NA19238.chrom22.SLX.maq.SRP000032.2009_07.bam
#	NA12004.SLX.maq.SRP000031.2009_09.unmapped.bam
#	NA10851.SOLID.SRP000031.2009_08.unmapped.bam  - does not have algorithm; not worry  since it is not going to be in BAM collection
	
#	NA12878.chrom1.LS454.ssaha.CEU.high_coverage.20091216.bam  - newest convention that includes the analysis group concept 
		
#	In index files, platform formats are:	
#	ABI_SOLID
#	ILLUMINA
#	LS454

sub CHECK_AND_PARSE_FILE_NAME {
	my ($file_path) = @_;
	
	my $file_name = basename($file_path);
	my @segments = split (/\./, $file_name);
	
	#print "file name is $file_name\n" if ($verbose);
	my $chr = 0;  ## this is to take care of BAMs that have all chromosomes merged
	my ($sample, $platform, $algorithm, $project, $date, $pop, $analysis_grp);
	if ( $file_path =~ /\/pilot_data\/data\// ) {  
				
		if ( $file_name =~ /chr/i  && @segments == 7 ) {
			($sample, $chr, $platform, $algorithm, $project, $date) = split(/\./, $file_name);
		}
		elsif ( $file_name =~ /unmapped/  &&  @segments == 7) {	
			($sample, $platform, $algorithm, $project, $date) = split(/\./, $file_name);
		}
		elsif (  $file_name =~ /unmapped/ &&  @segments == 6 ) {
			($sample, $platform, $project, $date) = split(/\./, $file_name);
		}	
		elsif ( @segments == 6 ) {
			($sample, $platform, $algorithm, $project, $date) = split(/\./, $file_name);
		}	
	}	
	else {  
		if (@segments == 7 && $file_name =~ /SRP/i) { ## main project naming convention, first version
			($sample, $chr, $platform, $algorithm, $project, $date) = split(/\./, $file_name);
		}
		elsif (@segments == 7) { ## main project convention, all chromosomes merged, NCBI bams
			($sample, $platform, $algorithm, $pop, $analysis_grp, $date) = split(/\./, $file_name);
		}	
		elsif (@segments == 8) { ## main project naming convention, latest version
			unless ($file_name =~ /SRP/i)  {
				($sample, $chr, $platform, $algorithm, $pop, $analysis_grp, $date) = split(/\./, $file_name);
			}
		}
		else {	
			throw("BAM file name $file_name does not fit naming convention\n");
		}	
	}
	
	### FIXME: need to do some controled vacabulary checks
	unless ( ($sample =~ /^NA/ || $sample=~/^HG/ || $sample =~ /^GM/) && (!$project || $project =~ /^SRP/) && $date =~ /2009|2010|2011|2012|2013/ ) {
		if ($project) {
			print "sample is $sample, project is $project, date is $date\n";
		}
		else {
			print "sample is $sample, date is $date\n";
		}	 	 
		throw("BAM file name does not fit naming convention\n");	
	}
	
	return ($sample, $platform, $algorithm, $project, $analysis_grp, $chr, $date, $pop);	
}		 

sub check_bas {
	
	my ($bas_path) = @_;
	
	my $flag = 0;
	
	open (BAS, "<", $bas_path) || throw("Cannot open bas file $bas_path\n");

	my $bas_basename = basename($bas_path);
	$bas_basename =~ s/\.bam\.bas//;
		
	while (<BAS>) {
		next if ($_ =~ /filename/);
		
		my @tmp = split (/\t/, $_);
		my $bam_file_name = $tmp[0];
		if ($bam_file_name ne $bas_basename) {
			$flag = 1;
		}	
			
	}
	close BAS;
	return $flag;	 
}


sub get_bam_run_id_info{

#BCM SOLID bams has RG in PU field of bam header. Must catch.
#@RG     ID:1    PL:SOLiD        PU:SRR097880    LB:ANG_TG.HG00142-1_1sA SM:HG00142      CN:BCM
#@RG     ID:2    PL:SOLiD        PU:SRR097881    LB:ANG_TG.HG00142-1_1sA SM:HG00142      CN:BCM
#
# Extract RG head info and change when in sub  correct_bas_file_convention

  my $bam  = shift;

  if ( !-e $bam){
    throw "$bam does not exist\n";
  }


  my $key = "-";
  my %rg_hash;

  my @rg_info =  `samtools view $bam -H | grep "^\@RG"`;


  die "No read group info found in $bam\n" if ( scalar  @rg_info < 1);
#  print "Found RG info for ",  scalar  @rg_info ," runs\n";

  foreach my $line (@rg_info) {
 #   print $line;
    my @aa = split /\s+/,$line;
    shift @aa;

    my $key = "-";
    foreach my $x (@aa) {
      my ($p1, $p2) = split /:/,$x;
      ($key = $p2) if ( $p1 eq "ID");
    }
    die "No ID in RG row \n" if ( $key eq "-");


    foreach my $x (@aa) {
      my ($p1, $p2) = split /:/,$x;
      next if ($p2 eq $key);
      $rg_hash{$key}{$p1} =  $p2;
    }
  }
 # print Dumper %rg_hash;

  return (\%rg_hash);
}
