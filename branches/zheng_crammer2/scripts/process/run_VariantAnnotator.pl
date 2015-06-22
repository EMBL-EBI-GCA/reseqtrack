#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use Getopt::Long;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use Time::Local;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::BamUtils;
use File::Path;
use File::Basename;
use ReseqTrack::Tools::RunVariantCall::RunVariantAnnotator;

my %input;

my $start_time = time();

&GetOptions( 
	\%input, 
	  	
	'dbhost=s',   
	'dbname=s',      
	'dbuser=s',
	'dbpass=s', 
	'dbport=s',  
  
	'program=s',
  	'reference:s',
  	'parameters=s%',
  	
	'input_bams:s@',
	'input_vcf=s',
	
	'anno:s@',
	'anno_group:s@',
	'dbsnp:s',

  	'chrom:s',
  	'region:s',
  	  	  	  	
  	'output_dir=s',
  	'output_file_type:s',

  	'store!',
 	'update!',
 	'host:s',
  	'save_files_from_deletion!',
  	'help!',
);	
	

if ( defined $input{cfg_file} ) {
  get_params( $input{cfg_file}, \%input );
}

### Set Default
$input{update} = 0 if (!$input{update});

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd` ;
	chomp $input{output_dir};
}

$input{output_dir} =~ s/\/$//;
mkpath($input{output_dir}) unless (-e $input{output_dir});

$input{save_files_from_deletion} 	= 0 if (!$input{save_files_from_deletion});	


=head
$input{host} 						= '1000genomes.ebi.ac.uk' if (!$input{host});
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $input{dbhost},
  -user => $input{dbuser},
  -port => $input{dbport},
  -dbname => $input{dbname},
  -pass => $input{dbpass},
);
    
my $fa = $db->get_FileAdaptor;
my $hist_a = $db->get_HistoryAdaptor;
my $ha = $db->get_HostAdaptor;

my $host = $ha->fetch_by_name($input{host});
if(!$host){
  $host = ReseqTrack::Host->new
      (
       -name => $input{host},
       -remote => $input{remote},
      );
}

$db->dbc->disconnect_when_inactive(1);
=cut
    		
my $object = ReseqTrack::Tools::RunVariantCall::RunVariantAnnotator->new(
	-program					=> $input{program},
	-reference 					=> $input{reference},
	-input_files				=> $input{input_bams},
	-input_vcf					=> $input{input_vcf},
	-anno 						=> $input{anno},
	-anno_group					=> $input{anno_group},
	-working_dir				=> $input{output_dir},
	-save_files_from_deletion	=> $input{save_files_from_deletion},
	-dbsnp						=> $input{dbsnp},
	-parameters					=> $input{parameters},
	-chrom						=> $input{chrom},
	-region						=> $input{region},
);

$object->run;

my $annotated_vcfs = $object->output_files;

foreach my $annotated_vcf ( @$annotated_vcfs  ) {  
	if ( -e $annotated_vcf ) {
		print "Annotated VCF file is $annotated_vcf\n";
	}
	else {
		print "Annotated VCF file $annotated_vcf doesn't exist\n";
	}
}	

my $end_time = time();
my $length = ($end_time - $start_time)/60; 
print "Job took $length minutes\n";

=pod

=head1 NAME

	~/ReseqTrack/scripts/process/run_VariantAnnotator.pl
	
=head1 SYNOPSIS

	This is a script to run Variant Annotator to add annotations to variant calls in a given VCF file.  
	Please see GATK website for detailed descriptions about variantAnnotator
	http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_annotator_VariantAnnotator.html
	
=head1 Arguments:

 Required arguments:
    
	reference,			Reference genome; 
	input_bams,			bam files that have been used to make the varaint calls
	input_vcf, 			input VCF file
	anno_group,			One or more classes/groups of annotations to apply to variant calls. Use the "-parameters ls=" argument to view available groups.
	anno,				One or more specific annotations to apply to variant calls. Use the "-parameters ls=" argument to view available annotations.

 Optional arguments:
 
 	parameters,			other parameters to pass on to VariantAnnotator, format is "-parameters parameter_name=parameter_value", can repeat for multiple entries 
 	dbsnp,
 	output_dir,
 	chrom,
 	region,
 	save_files_from_deletion

=head1	OUTPUT

	The output file of this script is a vcf file written to the output_dir; the file name is the input VCF file name with a suffix ".annotated.vcf" 

=head1 EXAMPLE COMMAND LINE	

perl $ZHENG_RT/scripts/process/run_VariantAnnotator.pl \
-reference /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/human_g1k_v37.fasta \
-dbsnp /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/dbsnp_135.b37.vcf \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_PUR_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_PUR_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_IBS_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CEU_LS454_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CHS_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_ASW_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_IBS_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CLM_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_MXL_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_YRI_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CHB_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_FIN_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_MXL_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_GBR_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CLM_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_LWK_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_YRI_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_JPT_SOLID_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_GBR_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_LWK_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_JPT_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CHB_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CHS_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_ASW_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_FIN_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_TSI_ILLUMINA_chrom20.bam \
-input_bams /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CEU_ILLUMINA_chrom20.bam \
-output_dir  /nfs/1000g-work/G1K/work/zheng/snp_calling/samtools \
-anno_group StandardAnnotation \
-input_vcf /nfs/1000g-work/G1K/work/zheng/snp_calling/samtools/all_chr20.new.q_gt_3.copy.vcf \
-chrom 20 \
-region 1000000-2000000 \




		