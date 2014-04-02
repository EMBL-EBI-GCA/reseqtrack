#!/usr/bin/env perl

use strict;
use warnings;
use ReseqTrack::Tools::Exception;
use Getopt::Long;
use Time::Local;
use File::Path;
use File::Basename;
use ReseqTrack::Tools::GATKTools::VariantEval;

my $start_time = time();

my %input;
my %hash;
my %comp_type_hash;
my %hash_to_print;
my $out_h;

&GetOptions( 
	\%input, 
	'program=s',
  	'reference:s',	
	'comp_vcfs:s@',
	'evalmodules:s@',
	'eval_vcf=s',
	'dbsnp:s',  	  	
  	'output_dir=s',
  	'help!',
);	

if ($input{help}) {
	 exec('perldoc', $0);
}

## Set GenotypeConcordance as one of the default evalModules
if ( $input{evalmodules} ) {
	my @tmp_array = qw(GenotypeConcordance);
	foreach my $mod (@{ref($input{evalmodules}) eq 'ARRAY' ? $input{evalmodules} : [$input{evalmodules}]}) {
    	push @tmp_array, $mod;  
    }
	$input{evalmodules} = \@tmp_array;
}
else {
	$input{evalmodules} =  qw(GenotypeConcordance);
}	

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd` ;
	chomp $input{output_dir};
}

$input{output_dir} =~ s/\/$//;
mkpath($input{output_dir}) unless (-e $input{output_dir});

my $outf = $input{output_dir} . "/" . basename($input{eval_vcf}) . ".summary"; 
open ($out_h, ">", $outf) || throw("Cannot open output file $outf");
  		
my $object = ReseqTrack::Tools::GATKTools::VariantEval->new(
	-program					=> $input{program},
	-reference 					=> $input{reference},
	-input_files				=> $input{eval_vcf},
	-comps						=> $input{comp_vcfs},
	-working_dir				=> $input{output_dir},
	-dbsnp						=> $input{dbsnp},
	-evalModules				=> $input{evalmodules}
);

$object->run;

my $report = $object->output_files->[0]; #As test:  my $report = "/nfs/1000g-work/G1K/work/zheng/vcf_report/merged.no_select.compOmni.eval.gatkreport";

throw("gatk report $report doesn't exist") unless ( -e $report );

$/ = "\n\n";

open(my $in, "<", $report) || throw("Cannot open GATK VCF report $report");

## Store the GATK report in a hash
while (<$in>) {
	chomp;
	my $block = $_;
	my @blk_data = split(/\n/, $block);	
	
	my @field_names;
	foreach my $line ( @blk_data ) {
		next if ($line =~ /^#/);
		my @line_data = split (/\s+/, $line);

		if ($line =~ /JexlExpression/) {
			@field_names = @line_data;
		}
		else {
			my $line_id = join("_", @line_data[0..4]);
			foreach my $index ( 5..$#line_data ) {
				$hash{$line_id}{$field_names[$index]} = $line_data[$index];
				#print $index . "\t" . $field_names[$index] . "\t" . $line_data[$index] . "\n";
			}	
		}	
	}									
}	

## Figure out how many types of comparison exist (none, dbSNP, OMNI, HapMap etc.) 
foreach my $key (keys %hash) {
	my ($analysis_type, $compRod, $evalRod, $JexlExp, $novelty) = split(/_/, $key);	
	$comp_type_hash{$compRod} = 1;
}	

## Pick out useful information from the GATK report, store them in a hash to print 
foreach my $key (keys %hash) {
	my ($analysis_type, $compRod, $evalRod, $JexlExp, $novelty) = split(/_/, $key);	
	if ( $analysis_type =~ /CountVariants/ || $analysis_type =~ /VariantSummary/ || $analysis_type =~ /CompOverlap/) {	
		foreach my $key2 ( keys %{$hash{$key}} ) {
			if ( 	( keys %comp_type_hash == 1 ) ||
					( keys %comp_type_hash >= 2 && $compRod =~ /comp/ ) ) {	
						## if no comp file (comp or dbSNP) was provided, store all rows 
						## OR if only one comp file (either comp or dbSNP) was provided, store all rows  
						## OR if two or more than two comp files (comp and dbSNP) were provided, only store the row(s) with comp or comp2 in the hash, ignore dbSNP rows
				$hash_to_print{$compRod}{$novelty}{$key2} = $hash{$key}{$key2};
			}
		}
	}		
} 
	
my @table1_headers = (	"nSamples",
						"nSNPs",
						"TiTvRatio",
						"nIndels",
						"nSVs",
						"nHets",
						"nHomRef",
						"nHomVar",
						"nSingletons",
						);	
														
unless (keys %comp_type_hash == 1 && defined $comp_type_hash{"none"} ) {  	# print the very basic information if no comparison files was provided, 
																			# print more columns if one or more comparison files were provided 
	push @table1_headers, 	(
							"novelSites",
							"nVariantsAtComp",
							"compRate",
							"nConcordant",
							"concordantRate"
							);	
}						
						
print $out_h "Comp\tNovelty\t" . join("\t", @table1_headers) . "\n";
foreach my $compRod ( keys %hash_to_print ) {
	foreach my $nov ( ("all", "known", "novel") ) {	
		print $out_h "$compRod\t$nov\t";			 
		foreach my $f (@table1_headers) {
			print $out_h $hash_to_print{$compRod}{$nov}{$f} . "\t";
		}
		print $out_h "\n";
	}		
}
	 
my $end_time = time();
my $length = ($end_time - $start_time)/60; 
print "Job took $length minutes\n";

=pod

=head1 NAME

	~/ReseqTrack/scripts/process/dump_vcf_report.pl
	
=head1 SYNOPSIS

	This is a script to run GATK's VariantEval and parse the resulted report to generate a basic VCF report.   
	Please see GATK website for detailed descriptions about VariantEval
	
	http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_varianteval_VariantEval.html#--comp
	http://www.broadinstitute.org/gatk/guide/article?id=2361
	
	Note that the GATK we installed in house is an old version; the GATK report may be slightly different from what is described in the web site.
	
=head1 Arguments:

 Required arguments:
    
	reference,			Reference genome
	eval_vcf, 			input VCF file to be evaluated
	program,			Path to the GATK program directory (Path to GATK program directory may also be provided through enviroment variable GATK)

 Optional arguments:
 
 	comp_vcfs,			input comparison files that contain datasets you want to compare your evaluation file to. This tag can be used multiple times.
 	dbsnp,				dbSNP file, this is used to classify the evaluation calls into "known" and "novel".
 	evalModules,		additional evalModules one wish to use.  Use the --list option in GATK VariantEval to see all possible modules (GenotypeConcordance is 
						default here in addition to other default modules set by GATK) 
 	output_dir,			Where temparary file and final output file will be written to.
  	help,
  	
=head1	OUTPUT

	The output file of this script is a vcf summary report. The name is the input evaluation vcf file with "summary" as suffix. 
	Below is an example of the report:
	
Comp  Novelty nSamples        nSNPs   TiTvRatio       nIndels nSVs    nHets   nHomRef nHomVar nSingletons     novelSites      nVariantsAtComp compRate        nConcordant     concordantRate
dbsnp   all     2       55403   2.47    0       0       45556   0       35455   22046   1583    53820   97.14   53722   99.82
dbsnp   known   2       53820   2.56    0       0       43893   0       35289   20820   0       53820   100.00  53722   99.82
dbsnp   novel   2       1583    0.91    0       0       1663    0       166     1226    1583    0       0.00    0       0.00
	
	The followings are the column headers of the VCF report:
	
	Basic columns:
	
	Comp,				This can be "comp", "comp1", or "comp2"; when no comp_vcfs are provided, this can be "dbSNP"; it will be "none" if even dbSNP file is not provided
	Novelty,			"all", "known", "novel" based on dbSNP
	nSamples,			Number of samples in the evaluation VCF file
	nSNPs,				Number of SNPs in the evaluation VCF file
	TiTvRatio,			Ti/Tv ratio
	nIndels,			Number of Indels in the evaluation VCF file 
	nSVs,				Number of SVs in the evaluation VCF file
	nHets,				Number of heterozygous genotypes in the evaluation VCF file
	nHomRef,			Number of homozygou ref genotypes in the evaluation VCF file
	nHomVar,			Number of homozygous alternative allele genotypes in the evaluation VCF file
	nSingletons,		Number of sites with alternative allele in a single sample
							
	Additional columns when comparison file(s) are provided with -comp_vcfs:
	 
	novelSites,			Number of sites only exist in the evaluation VCF file
	nVariantsAtComp,	Number of sites found in both evaluation and comparison VCF files
	compRate,			Percentage of sites in the evaulation file that are also found in comparison file 
	nConcordant,		number of concordant sites (that is, for the sites that share the same locus as a variant in the comp track, those that have the same alternate allele)
	concordantRate		Percentage of concordant sites
						  

=head1 EXAMPLE COMMAND LINE	

perl $ZHENG_RT/scripts/process/dump_vcf_report.pl \
-program /nfs/1000g-work/G1K/work/bin/gatk/dist \
-reference /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/human_g1k_v37.fasta \
-dbsnp  /nfs/1000g-archive/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz  \
-eval_vcf /nfs/1000g-work/G1K/work/zheng/vcf_report/CRI4237.vcf.gz \
-comp_vcfs  /nfs/1000g-archive/vol1/ftp/technical/working/20131122_broad_omni/Omni25_genotypes_2141_samples.b37.v2.vcf.gz \
-comp_vcfs /nfs/1000g-archive/vol1/ftp/technical/working/20110322_hapmap3_grch37/hapmap3.3.genotypes.b37_fwd.vcf.gz \
-output_dir ~ \
-evalmodules GenotypeConcordance \

							