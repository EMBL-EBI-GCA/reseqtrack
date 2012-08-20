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
use ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator;
use ReseqTrack::Tools::RunVariantCall::RunApplyRecalibration;

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
	'input_files=s',
	  	
  	'parameters_VR=s%',
	'resources=s%',
	'use_annotation=s@',
	
	'parameters_AR=s%',
	
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
	
### Set Default
$input{update} = 0 if (!$input{update});

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd` ;
	chomp $input{output_dir};
}

$input{output_dir} =~ s/\/$//;
mkpath($input{output_dir}) unless (-e $input{output_dir});

$input{save_files_from_deletion} 	= 0 if (!$input{save_files_from_deletion});	


#$input{host} 						= '1000genomes.ebi.ac.uk' if (!$input{host});
#$input{output_file_type} 			= "VQSR_FILTERED_VCF" if ( !$input{output_file_type} );

#my $db = ReseqTrack::DBSQL::DBAdaptor->new(
#  -host => $input{dbhost},
#  -user => $input{dbuser},
#  -port => $input{dbport},
# -dbname => $input{dbname},
# -pass => $input{dbpass},
#);
    
#my $fa = $db->get_FileAdaptor;
#my $hist_a = $db->get_HistoryAdaptor;
#my $ha = $db->get_HostAdaptor;

#my $host = $ha->fetch_by_name($input{host});
#if(!$host){
#  $host = ReseqTrack::Host->new
#      (
#       -name => $input{host},
#       -remote => $input{remote},
#      );
#}

#$db->dbc->disconnect_when_inactive(1);
   		
my $object_VR = ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator->new(
	-program					=> $input{program},
	-reference 					=> $input{reference},
	-input_files				=> $input{input_files},
	-resources 					=> $input{resources},
	-working_dir				=> $input{output_dir},
	-save_files_from_deletion	=> $input{save_files_from_deletion},
	-use_annotation				=> $input{use_annotation},
	-parameters_VR				=> $input{parameters_VR},
);

$object_VR->run;

my $intermediate_files = $object_VR->created_files;

my $tranchesFile;
my $recalFile;
foreach my $intermediate_file ( @$intermediate_files ) {  
	if ( -e $intermediate_file && $intermediate_file =~ /tranches/i ) {
		$tranchesFile = $intermediate_file;
	}
	elsif ( -e $intermediate_file && $intermediate_file =~ /recal/i ) {
		$recalFile = $intermediate_file;
	}
	else {
		print "intermediate file $intermediate_file doesn't exist\n";
	}
}	

my $object_AR = ReseqTrack::Tools::RunVariantCall::RunApplyRecalibration->new(
	-program					=> $input{program},
	-reference 					=> $input{reference},
	-input_files				=> $input{input_files},
	-working_dir				=> $input{output_dir},
	-save_files_from_deletion	=> $input{save_files_from_deletion},
	-tranchesFile				=> $tranchesFile,
	-recalFile					=> $recalFile,
	-parameters_AR				=> $input{parameters_AR},
);

$object_AR->run;

my $calibrated_vcf = $object_AR->output_files;

foreach my $outfile ( @$calibrated_vcf ) {  
	if ( -e $outfile ) {
		print "VQSR calibrated output VCF file is $outfile\n";
	}
	else {
		print "VQSR calibrated output VCF file $outfile doesn't exist\n";
	}
}	

my $end_time = time();
my $length = ($end_time - $start_time)/60; 
print "Job took $length minutes\n";


=pod

=head1 NAME

	~/ReseqTrack/scripts/process/run_VQSR.pl
	
=head1 SYNOPSIS

	This is a script to run Variant Quality Score Recalibrator (VQSR) to refine/filter variant calls made from various calling programs.  
	Please see GATK website for detailed descriptions about VQSR
	http://gatkforums.broadinstitute.org/discussion/39/variant-quality-score-recalibrator
	
=head1 Arguments:

 Required arguments:
    
	reference,			Reference genome; 
	resources,			Training sets of high quality known variants; can be multiple sets such as hapmap, omni variants and dbSNP
						Should be provided in format of:
						resource_name="file features, file path"
						"File features" is to indicate whether the file is to be used as "known", "training" or "truth" set with a "prior" probability of being correct
						Example: 
						dbSNP="known=true,training=false,truth=false,prior=8.0 /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/dbsnp_135.b37.vcf"
	input_files,		Raw variant calls that need to be recalibrated
	use_annotation, 	The names of the annotations which should used for calculations (QD, HaplotypeScore, MQRankSum, ReadPosRankSum, MQ etc)
	parameters_VR,		other parameters to pass on to VariantRecalibrator, format is parameter_name=parameter_value 
	parameters_AR, 		other parameters to pass on to ApplyRecalibration, format is parameter_name=parameter_value 

 Optional arguments:
 
 	output_dir,
 	save_files_from_deletion

=head1	OUTPUT

	The output file of this script is a vcf file written to the output_dir; the file name is the input VCF file name with a suffix ".vqsr_filtered.vcf" 

=head1 EXAMPLE COMMAND LINE	

perl $ZHENG_RB_VQSR/scripts/process/run_VQSR.pl \
-reference /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/human_g1k_v37.fasta \
-resources dbSNP="known=true,training=false,truth=false,prior=8.0 /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/dbsnp_135.b37.vcf" \
-resources hapmap="known=false,training=true,truth=true,prior=15.0 /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/hapmap_3.3.b37.sites.vcf" \
-resources omni="known=false,training=true,truth=false,prior=12.0 /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/1000G_omni2.5.b37.sites.vcf" \
-input_files /nfs/1000g-work/G1K/work/zheng/snp_calling/gatk/gatk_all_chr20.vcf.gz \
-output_dir /nfs/1000g-work/G1K/work/zheng/snp_calling/vqsr/results \
-use_annotation QD \
-use_annotation HaplotypeScore \
-use_annotation MQRankSum \
-use_annotation ReadPosRankSum \
-use_annotation MQ \
-parameters_VR -maxGaussians=6 \
-parameters_VR mode=BOTH \
-save_files_from_deletion \
-parameters_AR -ts_filter_level=99.0 \
-parameters_AR mode=BOTH \