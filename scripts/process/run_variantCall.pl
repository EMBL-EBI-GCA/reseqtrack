#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::GeneralUtils;
use Getopt::Long;
use ReseqTrack::Tools::BamUtils;
use File::Basename;
use ReseqTrack::Tools::Loader::File;

#use ReseqTrack::Tools::RunVariantCallUtils qw(validate_input_string);

$| = 1;

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
	
	$collection_name,
    $collection_type,
    $file_name_pattern,
	$bam_list,

	$reference,	
	$module_name,
	$program_path,
    $chrom,
    $region,
    $algorithm,
    $parameters,
	
	@bams,
	$verbose,
	$store,
	$help,
   );
 
my $update = 0;
my $output_dir = `pwd`;
chomp $output_dir;
my $host_name = '1000genomes.ebi.ac.uk';
my $output_file_type = "VCF";
my $module_path = 'ReseqTrack::Tools::RunVariantCall';
my $keep_intermediate_file = 0;

&GetOptions( 
  'dbhost=s'      		=> \$dbhost,
  'dbname=s'      		=> \$dbname,
  'dbuser=s'      		=> \$dbuser,
  'dbpass=s'      		=> \$dbpass,
  'dbport=s'      		=> \$dbport,
  
  'collection_name=s' 	=> \$collection_name,
  'collection_type=s' 	=> \$collection_type,
  'file_name_pattern=s' => \$file_name_pattern,
  'bam_list=s'			=> \$bam_list,
  
  'program=s'			=> \$program_path,
  'output_dir=s'		=> \$output_dir,
  'reference:s' 		=> \$reference,
  'chrom:s'				=> \$chrom,
  'region:s'			=> \$region,
  'algorithm=s'			=> \$algorithm,
  'parameters|options=s'=> \$parameters,
  'output_file_type:s'	=>\$output_file_type,
  'store!' 				=> \$store,
  'update!' 			=> \$update,
  'host:s'				=> \$host_name,
  'keep_intermediate_file!'	=>\$keep_intermediate_file,
  'help!'				=> \$help,
);

if ($help) {
	help_info();
}

if ( !$algorithm || $algorithm !~ /samtools|umake|gatk/i ) {
	throw("Please provide snp calling algorithm you wish to use, can be samtools, gatk, or umake\n");
}

if ( !$reference) {
	warn("No reference genome is provided, will use the default reference genome");
}	

if(!($collection_name && $collection_type) && !$bam_list && !$file_name_pattern){
  throw("Please provide some input arguments either -collection_name and ".
        "-collection_type, -bam_list or -file_name_pattern");
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $fa = $db->get_FileAdaptor;

if($collection_type && $collection_name){
  my $collection =
    $db->get_CollectionAdaptor->fetch_by_name_and_type($collection_name,
                                                        $collection_type);
  foreach my $file (@{$collection->others}){
    next if ( check_file_sanity($file->name, $chrom) eq 'skip' );
    push(@bams, $file->name);
  }
}

if($bam_list){
  my $lines = get_lines_from_file($bam_list);
  foreach my $line(@$lines){
    next if ( check_file_sanity($line, $chrom) eq 'skip');
    push(@bams, $line);
  }
}

if($file_name_pattern){
  my $files = $db->get_FileAdaptor->fetch_all_like_name($file_name_pattern);
  foreach my $file(@$files){
    next if ( check_file_sanity($file->name, $chrom) eq 'skip' );
    push(@bams, $file->name);
  }
}

$db->dbc->disconnect_when_inactive(1);

print "Input BAMs are:\n";
print join ("\n", @bams), "\n";

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

my $varCalling_module = $module_path . '::'. $module_name;

my $constructor_hash = parameters_hash($parameters, $algorithm);
$constructor_hash->{-input_files} = \@bams;
$constructor_hash->{-program} = $program_path;
$constructor_hash->{-working_dir} = $output_dir;
$constructor_hash->{-reference} = $reference;
$constructor_hash->{-chrom} = $chrom if ($chrom);
$constructor_hash->{-region} = $region if ($region);
$constructor_hash->{-save_files_for_deletion} = $keep_intermediate_file;

my $object = construct_new_object($varCalling_module, $constructor_hash);

$object->run;

my $flt_vcf = $object->output_files;

$db->dbc->disconnect_when_inactive(0);

foreach my $fo ( @$flt_vcf ) {  
	if ( -e $fo ) {
				
		my $loader = ReseqTrack::Tools::Loader::File->new
		  (
		   -file => [$fo],
		   -do_md5 => 1,
		   -hostname => $host_name,
		   -db => $db,
		   -assign_types => 0,
		   -check_types => 0,
		   -type => $output_file_type,
		   -update_existing => $update,
		  );
		
		if($store){
		  $loader->process_input();
		  $loader->create_objects();
		  $loader->sanity_check_objects();
		  my $file_objects = $loader->load_objects();
		}
		
		print "*********************************************************************************************************************************\n";
		print "**** Output filtered VCF file path is $fo\n";
		print "**** Output filtered VCF file $fo has been stored in the database\n" if ($store);
		print "*********************************************************************************************************************************\n";  
	}	
	else {
		throw("Output file " . $fo . " does not exist");
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
#          print "key is $key, value is $value\n";
#         validate_input_string($key, $algorithm);		 
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

sub check_file_sanity{
  my ($file, $chrom) = @_;
  my $flag = 'ok';
 
  throw($file." does not exist on the current filesystem") unless(-e $file);
 
  $flag = 'skip' if ($file !~ /\.bam$/);
 
  if ($chrom) {
      my $chr = "chrom" . $chrom;
      $flag = 'skip' if ( $file !~ /$chr/ );
  }     
  return $flag;
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
			 
	-collection_name		Name of a BAM collection
	-collection_type		Type of BAM collection
	-file_name_pattern		Any text string that exist in file names
	-bam_list,			Instead of providing sample names, a list of paths to BAM files can be used as input
	-chrom,				optional for samtools and gatk but mandatory for umake
	-region,			Format: 1000000-2000000
	-reference,			Reference genome; default is ncbi build 37  
	-output_dir,			This is where all temporary files and output files will be written to
	-store,				To store the output VCF file in the database
	-update,			To replace an output VCF file in the database
	-keep_intermediate_file		Do not delete intermediate files
	
	-parameters,			Key and value pairs to pass options to different variant calling algorithms. Options are different for each algorithm:
	-options,			same as -parameters

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
-ref /nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta \
-bam_list /nfs/1000g-work/G1K/work/zheng/snp_calling/bam_list \
-parameters 'mpileup=>-ug,bcfview=>-bvcg,vcfutils=>-D 100' \
-store \
-update	&

perl ~/ReseqTrack/scripts/process/run_variantCall.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-algorithm samtools \
-ref /nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta \
-file_name_pattern NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.2011b \
-chrom 11 \
-parameters 'mpileup=>-ug,bcfview=>-bvcg,vcfutils=>-D 100' &
	
An example that is quick to run:	
perl ~/ReseqTrack/scripts/process/run_variantCall.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-algorithm samtools \
-ref /nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta \
-collection_name NA06985.ILLUMINA.bwa.low_coverage \
-collection_type TEST_BAM \
-parameters 'mpileup=>-ug,bcfview=>-bvcg,vcfutils=>-D 100' \
-chrom 10 \
-region 1000000-2000000 &	

=head1 GATK EXAMPLES

perl ~/ReseqTrack/scripts/process/run_variantCall.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-algorithm gatk \
-ref /nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta \
-bam_list /nfs/1000g-work/G1K/work/zheng/snp_calling/bam_list \
-chrom 10 -region 10000000-15000000 \
-parameters 'dcov=>10,stand_emit_conf=>10.0,stand_call_conf=>50.0,glm=>BOTH,dbSNP=>/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/dbsnp132_20101103.vcf.gz' &


=head1 UMAKE EXAMPLES

perl ~/ReseqTrack/scripts/process/run_variantCall.pl \
-dbhost mysql-g1kdcc-public -dbname g1k_archive_staging_track -dbuser g1krw -dbpass thousandgenomes -dbport 4197 \
-algorithm umake \
-bam_list /nfs/1000g-work/G1K/work/zheng/snp_calling/bam_list \
-parameters 'dbSNP=>/nfs/1000g-work/G1K/work/bin/umake-resources/dbSNP/dbsnp_129_b37.rod,hm3_prefix=>/nfs/1000g-work/G1K/work/bin/umake-resources/HapMap3/hapmap3_r3_b37_fwd.consensus.qc.poly,indel_prefix=>/nfs/1000g-work/G1K/work/bin/umake-resources/indels/1kg.pilot_release.merged.indels.sites.hg19' \
-reference /nfs/1000g-work/G1K/work/bin/umake-resources/ref/human.g1k.v37.fa \
-chrom 10 -region 1-10000000 \
-output_dir /nfs/1000g-work/G1K/work/zheng/snp_calling/umake/1103 & 
 
 
