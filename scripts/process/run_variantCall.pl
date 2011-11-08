#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use Getopt::Long;
use ReseqTrack::Tools::BamUtils;
use File::Basename;

use ReseqTrack::Tools::RunVariantCallUtils qw(validate_input_string);

$| = 1;

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
	
	$sample_string,
	$bam_type,
	$analysis_grp,
	$seq_platform,

	$bam_list,
	$program_path,
    $chrom,
    $region,
    $algorithm,
    $parameters,
	
	$verbose,
	$store,
	$help
   );
 
my $reference = "/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta";    
my $update = 0;
my $output_dir = `pwd`;
chomp $output_dir;
my $host_name = '1000genomes.ebi.ac.uk';

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  
  'samples=s'		=> \$sample_string, # $samples (samples separated by , or space -- NA12345,NA12456)
  'bam_type=s'		=> \$bam_type,
  'analysis=s'		=> \$analysis_grp,
  'platform=s'		=> \$seq_platform,
  
  'bam_list=s'		=> \$bam_list,
  'program=s'		=> \$program_path,
  'output_dir=s'	=> \$output_dir,
  'reference:s' 	=> \$reference,
  'chrom:s'			=> \$chrom,
  'region:s'		=> \$region,
  'algorithm=s'		=> \$algorithm,
  'parameters:s'	=> \$parameters,
  'store!' 			=> \$store,
  'update!' 		=> \$update,
  'host:s'			=> \$host_name,
  'help!'			=> \$help,
);

if ($help) {
	help_info();
}

if ( !$algorithm || $algorithm !~ /samtools|umake|gatk/i ) {
	throw("Please provide snp calling algorithm you wish to use, can be samtools, gatk, or umake\n");
}

if ( $algorithm =~ /umake/i && ! $chrom) {
	throw("When run umake, please specify a chromosome");
}	

if (!$sample_string && !$bam_list) {
	throw("Please provide sample name(s) and associated BAM type, analysis type and seq platform type; or provide a list of bam files.Examples:
-bam_list /nfs/1000g-work/G1K/work/zheng/snp_calling/bam_list 
or
-samples NA19240.chrom22 -bam_type BAM -analysis high_coverage -platform ILLUMINA \n");
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $fa = $db->get_FileAdaptor;

my 	@bams;
if ( $sample_string ) {

	my @samples;
	if ($sample_string =~ /,/ ) {
		@samples = split(/,/, $sample_string);
	}
	elsif ($sample_string =~ /\s/) {
		@samples = split(/\s+/, $sample_string);
	}		
	else {		## when input is a single sample
		@samples = ($sample_string);
	}
	
	my @files;
	foreach my $sample ( @samples ) {	 
		@files = @{$fa->fetch_all_like_name($sample)};
		foreach my $file ( @files ) {
			next if ( $file->type ne $bam_type );
			my $bam_base = basename($file->name);
			my ($collection_name, $sam, $platform, $algo, $project, $analysis, $chr, $date) 
				= get_collection_name_from_file_name($bam_base);
			#print "analysis in file name is $analysis\nanalysis in input is $analysis_grp\n";	
			next if ($analysis !~ $analysis_grp) ;	
			next if ($platform !~ $seq_platform);
			push @bams, $file->name;
		}
	}		 
}

if($bam_list){
  my $lines = get_lines_from_file($bam_list);
  foreach my $path (@$lines){
    throw("Can't process ".$path." as it doesn't exist") unless(-e $path);
    push(@bams, $path);
  }
}

print "Input BAMs are:\n";
print join ("\n", @bams), "\n";

my $module_name;
if ($algorithm =~ /samtools/i) {
	$module_name = "CallBySamtools";
}	
elsif ( $algorithm =~ /gatk/i || $algorithm =~ /unified_genotyper/i ) {
	$module_name = "CallByGATK";
}
elsif ( $algorithm =~ /umake/i ) {
	$module_name = "CallByUmake";
}	
else {
	throw("Please use one of the 3 algorithms to make variant calls - samtools, gatk or umake\n");
}

my $varCalling_module = 'ReseqTrack::Tools::RunVariantCall::'.$module_name;

my $constructor_hash = parameters_hash($parameters, $algorithm);
$constructor_hash->{-input_files} = \@bams;
$constructor_hash->{-program} = $program_path;
$constructor_hash->{-working_dir} = $output_dir;
$constructor_hash->{-reference} = $reference;
$constructor_hash->{-chrom} = $chrom if ($chrom);
$constructor_hash->{-region} = $region if ($region);

my $object = construct_new_object($varCalling_module, $constructor_hash);

$object->run;

my $flt_vcf = $object->output_flt_vcf;
my $host = get_host_object($host_name, $db);
my $flt_vcf_obj = create_objects_from_path_list($flt_vcf, 'VCF', $host); 

foreach my $fo ( @$flt_vcf_obj ) {  
	if ( -e $fo->name ) {
		my $md5 = run_md5($fo->name) if ($store);
		$fo->md5($md5);
		$fa->store($fo, $update) if ($store);
		print "*****************************************************************\n";
		print "Output filtered VCF file path is " . $fo->name . "\n";
		print "Output filtered VCF file " . $fo->name . "has been stored in the database\n" if ($store);
		print "*****************************************************************\n";  
	}	
	else {
		throw("Output file " . $fo->name . " does not exist");
	}	
}


##############
#### SUBS ####
##############

sub parameters_hash{
  my ($string, $algorithm2) = @_;

  my %parameters_hash;

  if ($string) {
      
    my @pairs;
    if ($string !~ /,/ && $string =~ /=>/ ) {
		push @pairs, $string;
    }       
    elsif($string =~  /,/ && $string =~ /=>/){
		@pairs = split (/,/, $string);
    }          
    else {
          throw("Please provide both parameter name and value linked by =>");
	}
	
	foreach my $pair(@pairs){
		my ($key, $value) = split (/=>/, $pair);
          $key   =~ s/^\s+//g;
          $key   =~ s/\s+$//g;
          $value =~ s/^\s+//g;
          $value =~ s/\s+$//g;
          $key = "-" . $key;
          $parameters_hash{$key} = $value;
          print "key is $key, value is $value\n";
          validate_input_string($key, $algorithm);		 
    }#end of foreach pairs
  }# end of if ($string)
  
  return \%parameters_hash;
}  		 		
  		
sub construct_new_object {
  my ($module, $args) = @_;
  my $file = "$module.pm";
  $file =~ s{::}{/}g;
  eval {
    require "$file"; # raises exception upon error, checks @INC
  };
  if($@){
    throw("ReseqTrack::Tools::EventPipeline::setup_batch_submission_system ".
          "Can't find $file [$@]");
  }
  my %constructor_args = %$args;
  my $object = $module->new(%constructor_args);
  return $object;
}
			    
sub help_info {

 exec('perldoc', $0);

}


#############################################################################################
=pod

=head1 NAME

	~/ReseqTrack/scripts/process/run_variantCall.pl
	
=head2 Required arguments:

	-dbhost, 			the name of the mysql-host
	-dbname, 			the name of the mysql database
	-dbuser, 			the name of the mysql user
	-dbpass, 			the database password if appropriate
	-dbport, 			the port the mysql instance is running on, this defaults to 4197 the standard
						port for mysql-g1kdcc.ebi.ac.uk
    
	-algorithm,			can be one of the three: samtools, gatk, umake
	-chrom,				for umake, chr is mandatory; otherwise the makefile would have nothing. For samtools and gatk, it is optional
	
=head2 Optional arguments:		

	-samples,			one or more sample ids separated by , or space
	-bam_type,			type of BAM
	-analysis,			can be one of the four: low_coverage, high_coverage, exome, exon_targetted
	-paltform,			can be one of the three: ILLUMINA, SOLID, 454 (check exact words)			 

	-bam_list,			Instead of providing sample names, a list of paths to BAM files can be used as input

	-chrom,				optional for samtools and gatk but mandatory for umake
	-region,			Format: 1000000-2000000
	-reference,			Reference genome; default is ncbi build 37  
	-output_dir,		This is where all temporary files and output files will be written to
	-store,				To store the output VCF file in the database
	-update,			To replace an output VCF file in the database
	
	-parameters,		Key and value pairs to pass options to different variant calling algorithms. Options are different for each algorithm:
	
		Samtools:
		
		parameter keys have to be one or more of the followings. Please see Samtools page for details -
		http://samtools.sourceforge.net/samtools.shtml#3
		
		mpileup=>, 
		bcfview=>, 
		vcfutils=>, 
		bcftools_path=>, 
		vcfutils_path=>
		
		For mpileup=>, values can be taken from:
		samtools mpileup [-EBug] [-C capQcoef] [-l list] [-M capMapQ] [-Q minBaseQ] [-q minMapQ]
		The [-r region] and [-f reference] have been coded as separate input entries so don't include them in the parameteres
		
		For bcfview=>, values can be taken from:
		bcftools view [-AbFGNQSucgv] [-D seqDict] [-l listLoci] [-i gapSNPratio] [-t mutRate] [-p varThres] [-P prior] [-1 nGroup1] [-d minFrac] [-U nPerm] [-X permThres] [-T trioType]  
		[-s listSample] and [region] have been coded separately so don't include them in the parameters
		
		For vcfUtils=>, values can be
		The -D option of varFilter controls the maximum read depth, which should be adjusted to about twice the average read depth. 
		
		gatk:
		
		Parameter keys have to be one or more of the followings. Please see GATK website for detailed descriptions - 
		http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--output_mode
		dbSNP=> 
  		dcov=> 
  		computeSLOD=>
  		stand_emit_conf=> 
  		stand_call_conf=> 
  		debug_file=>
  		glm=>
  		genotyping_mode=> 
  		pcr_error_rate=> 
  		p_nonref_model=> 
  		minIndelCnt=>
  		output_mode=>
  		min_mapping_quality_score=>
  		min_base_quality_score=>
  		metrics_file=>
  		max_deletion_fraction=>
  		indel_heterozygosity=>
  		heterozygosity=>
  		group=>
  		annotation=>
  		alleles=>
		
		umake:

		Parameters keys have to be one or more of the following. Details can be found in UMAKE website - 
		http://genome.sph.umich.edu/wiki/UMAKE#Configuration_File
				
		offset_off_target=>
		target_bed=>  
  		dbSNP=> 
  		indel_prefix=> 
  		hm3_prefix=> 
  		FILTER_MAX_SAMPLE_DP=> 
  		FILTER_MIN_SAMPLE_DP=>
		
		
=head1 SYNOPSIS
	
	This is a generic script to call variants from BAM files. Algorthms that can be used for making varaint call include: 
	samtools (http://samtools.sourceforge.net/mpileup.shtml)   
	gatk (http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper)
	umake (http://genome.sph.umich.edu/wiki/UMAKE)
	
	Input BAMs are defined either by a list of file paths, or sample names, BAM type, analysis group, and sequence platform, 
	which can be used to query the database to get a list of BAM file path.  Output is variants discovered in the format of VCF file. 
	
=head1	OUTPUT
	
	Output files from samtools and gatk can be found in the user-specified output_dir; if an output_dir is not specified, the current 
	directory is the output_dir.
	
	Output file name contains one of the samples, algorthm name, and chromosomal region if specified. Below is a few examples of VCF files:
	NA06985_and_0_others.samtools.flt.vcf
	NA06985_and_3_others.chr10_10000000-20000000.gatk.vcf
	
	Output for UMAKE is defined by UMAKE itself and is more complex. They are organized in different directories specified in section "RELATIVE DIRECTORY UNDER OUT_DIR" 
	of the configuration file and separated into different chromosomes. Please browse the output dir for useful outputs.
	For --snpcall, the VCF file can be found under the 'vcfs/chrx' folder.
	
	The script run_variantCall.pl will print the path to the most relevant output file when a run is done.  Output VCF files will be stored 
	in the database as type "VCF" if option '-store' is specified.
	
=head1 SAMTOOLS EXAMPLES

perl ~/ReseqTrack/scripts/process/run_variantCall.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-algorithm samtools \
-bam_list /nfs/1000g-work/G1K/work/zheng/snp_calling/bam_list \
-parameters 'mpileup=>-ug,bcfview=>-bvcg,vcfutils=>-D 100' \
-store \
-update	&

perl ~/ReseqTrack/scripts/process/run_variantCall.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-algorithm samtools \
-samples NA19240.chrom22 \
-bam_type BAM \
-analysis high_coverage \
-platform ILLUMINA \
-parameters 'mpileup=>-ug,bcfview=>-bvcg,vcfutils=>-D 100' &
	
An example that is quick to run:	
perl ~/ReseqTrack/scripts/process/run_variantCall.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-algorithm samtools \
-samples NA06985.chrom10.ILLUMINA.bwa.CEU.low_coverage.2011b \
-bam_type TEST_BAM \
-analysis low_coverage \
-platform ILLUMINA \
-parameters 'mpileup=>-ug,bcfview=>-bvcg,vcfutils=>-D 100' \
-chrom 10 \
-region 1000000-2000000 &	
 	
=head1 GATK EXAMPLES

perl ~/ReseqTrack/scripts/process/run_variantCall.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-algorithm gatk \
-samples NA06985.chrom10.ILLUMINA.bwa.CEU.low_coverage.2011b \
-bam_type TEST_BAM \
-analysis low_coverage \
-platform ILLUMINA \
-chrom 10 -region 10000000-15000000 \
-parameters 'dcov=>10,stand_emit_conf=>10.0,stand_call_conf=>50.0,glm=>BOTH,dbSNP=>/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/dbsnp132_20101103.vcf.gz' &


=head1 UMAKE EXAMPLES
 
perl ~/ReseqTrack/scripts/process/run_variantCall.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-algorithm umake \
-samples NA12340.chrom20.ILLUMINA.bwa.CEU.low_coverage.20101123 \
-bam_type BAM \
-analysis low_coverage \
-platform ILLUMINA \
-parameters 'dbSNP=>/nfs/1000g-work/G1K/work/bin/umake-resources/dbSNP/dbsnp_129_b37.rod,hm3_prefix=>/nfs/1000g-work/G1K/work/bin/umake-resources/HapMap3/hapmap3_r3_b37_fwd.consensus.qc.poly,indel_prefix=>/nfs/1000g-work/G1K/work/bin/umake-resources/indels/1kg.pilot_release.merged.indels.sites.hg19' \
-reference /nfs/1000g-work/G1K/work/bin/umake-resources/ref/human.g1k.v37.fa \
-chrom 20 -region 1-10000000 \
-output_dir /nfs/1000g-work/G1K/work/zheng/snp_calling/umake/1103 & 
 
 