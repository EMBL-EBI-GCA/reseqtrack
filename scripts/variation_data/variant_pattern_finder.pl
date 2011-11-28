#!/sw/arch/bin/perl -w

=head1 LICENSE

	Copyright (c) 1999-2011 The European Bioinformatics Institute and
	Genome Research Limited.  All rights reserved.

	This software is distributed under a modified Apache license.
	For license details, please see

	http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

	Please email comments or questions to the 1000 genomes project information 
	list at <info@1000genomes.org>

=head1	AUTHOR

  	Holly Zheng Bradley (zheng@ebi.ac.uk)

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Vcf;

my (	$file,
		$region,					
		$sample_panel,
		%sample_pop,
		%local_snps,				# snps in a region
		%snp_impact,				# snps with non-synonymous or splice-site as consequences
		%snp_ref,					# key is snp identifier (pos + rs), value is REF allele
		%snps,						# key is the combo identifier for snp; value is the chr position
		%sample_gt,					# this stores any samples in a region and snp alleles as array, sorted by snp position
		%local_samples,				# this stores all samples with alt alleles in any snps in the region
		%unique_allele_strings,
		$help,
);

my $output_dir = ".";
my $print_all = 0;
my $expanded_view = 0;
my $verbose = 0;
	
&GetOptions( 
	'vcf=s'			=>\$file,
	'region=s'		=>\$region,
	'sample=s'		=>\$sample_panel,
	'output_dir=s'	=>\$output_dir,
	'expanded_view!'	=>\$expanded_view,
	'print_all!'	=>\$print_all,
	'help!'			=>\$help,
	'verbose!'		=>\$verbose,
);

if ($help) {
	 exec('perldoc', $0);
}

die("please provide vcf file, chromosomal region, and sample panel file") if (!$file || !$region || !$sample_panel);

my $outfile =  "$output_dir/chr" . $region . ".csv";
$outfile =~ s/:/_/;
open(OUT, ">", $outfile) || die("Cannot open output file $outfile");

my $header = "$output_dir/header";
open(HEADER, ">", $header) || die("Cannot open output header file $header");

my $total_sample_cnt = parse_sample_file($sample_panel);

my $vcf	= Vcf->new(file=>$file, region=>$region, print_header=>1); #print_header=>1 allows print sample name rather than column index
$vcf->parse_header();

while (my $x=$vcf->next_data_hash()) {
        #print join ("\t", keys (%$x)) , "\n"; # to get all possible hash keys:
        #FORMAT  QUAL    ID      CHROM   INFO    FILTER  gtypes  REF     ALT     POS 
        
        my $position = $x->{CHROM}.":".$x->{POS};
        print  "SNP: $position " . $x->{ID} . " REF IS: " . $x->{REF} . " ALT IS: " . $x->{ALT}->[0] . "\n" if ($verbose);

		my $identifier = $position . "_" . $x->{ID};
		
		my $func_flag = populate_snp_impact_hash($x, $identifier);
		
		next if ($func_flag == 0 && !$print_all);
		
        for my $individual (keys %{$$x{gtypes}}) {
            my ($al1,$sep,$al2) = $vcf->parse_alleles($x,$individual);
            
            #print "raw data: $individual: $al1$sep$al2\n";
            next if ($al1 eq $x->{REF} && $al2 eq $x->{REF});   
            next if ($al1 eq "." && $al2 eq "."  && !$expanded_view);      
                                                           
            my $pos_gt = "$position " . "$al1$sep$al2";        
            
            $local_snps{$identifier}{$individual} = "$al1$sep$al2";
            print "$individual: $al1$sep$al2\n" if ($verbose);
            $snp_ref{$identifier} = $x->{REF};

        }
}

#### Populate %snps, memorize all SNPs in the region, then sort them by position	####
#### Populate %local_samples, memorize all samples in the region					####
foreach my $snp1 (keys %local_snps) {

	my ($chr, $pos_rs) = split(/:/, $snp1);
	my ($pos, $rs) = split(/_/, $pos_rs);
	
	$snps{$snp1} = $pos;
	foreach my $sample1 ( keys %{$local_snps{$snp1}} ) {	
		$local_samples{$sample1} = 1;
	}	 		
}
my @sorted_snps = sort{$snps{$a}<=>$snps{$b}} keys %snps;  # sort hash keys based on hash values

#### Populate %sample_gt, memorize each sample's genotype ####
foreach my $snp ( @sorted_snps) {
	foreach my $sample ( keys %local_samples ) {	
		if ($local_snps{$snp}{$sample}) {
			push @{$sample_gt{$sample}}, $local_snps{$snp}{$sample};
		}
		else {
			push @{$sample_gt{$sample}}, "-" if ( !$expanded_view );
			push @{$sample_gt{$sample}}, "REF|REF" if ( $expanded_view );
		}
	}	
}		

#### Populate %unique_allele_strings, find out how many types of unique strings exist ####
foreach my $s (keys %sample_gt) {
	my $allele_string;
	foreach my $genotype ( @{$sample_gt{$s}} ) {
		$allele_string .= $genotype . "\t";
	}	
	$unique_allele_strings{$allele_string}{$s}=1;
}

#### Print SNP header ####
print HEADER "freq\t";
foreach my $sorted_snp ( @sorted_snps ) {
	print HEADER $sorted_snp . "[" . $snp_ref{$sorted_snp} . "]\t";
}
print HEADER "\tsamples\n";

#### Print SNP function header ####
print HEADER "freq\t";
foreach my $snp2 ( @sorted_snps ) {
	my @funcs = keys %{$snp_impact{$snp2}};
	if ( @funcs > 0 ) {
		print HEADER join(",",  @funcs) . "\t";
	}
	else {
		print HEADER "non_functional\t";
	}	
}
print HEADER "\tsamples\n";	
	
#### Print allele string for each sample combination ####
foreach my $allele_chain (keys %unique_allele_strings) {
	my @samples =  keys %{$unique_allele_strings{$allele_chain}};
	my %samples_to_print_by_pop;
	my %sam_cnt_by_pop;
	foreach my $samp ( @samples ) {
		push @{$samples_to_print_by_pop{$sample_pop{$samp}}}, $samp;
		if ( ! $sam_cnt_by_pop{$sample_pop{$samp}} ) {
			$sam_cnt_by_pop{$sample_pop{$samp}} = 1;
		} else {
			$sam_cnt_by_pop{$sample_pop{$samp}}++;
		}	
	}	
		
	my $frequency = sprintf ("%.2f", (@samples/$total_sample_cnt)*100 );
	print OUT "$frequency\t";
	print OUT $allele_chain;	
	
	my @sorted_sam_cnt_by_pop_keys = sort{$sam_cnt_by_pop{$b}<=>$sam_cnt_by_pop{$a}} keys %sam_cnt_by_pop;  # sort hash keys based on hash values
	foreach my $pop (@sorted_sam_cnt_by_pop_keys) {
		print OUT "$pop($sam_cnt_by_pop{$pop})\t";
		if ( $sam_cnt_by_pop{$pop} > 2) {
			my $number = $sam_cnt_by_pop{$pop} - 3;
			my @first_three_sams = ($samples_to_print_by_pop{$pop}[0], $samples_to_print_by_pop{$pop}[1], $samples_to_print_by_pop{$pop}[2]);
			print OUT join (", ", @first_three_sams) . " and $number others.\t";
		}
		else {	
			print OUT join(", ", @{$samples_to_print_by_pop{$pop}}) . "\t";
		}
	} 
	print OUT "\n";
}
close OUT;
close HEADER;

post_processing($outfile);

################
##### SUBS #####
################
sub parse_sample_file {
	my ($sample_panel) = @_;
	open(SAM, "<", $sample_panel) || die("Cannot open sample panel file $sample_panel"); ## FIXME: allow open a URL
	my %total_sample_cnt_hash;
	while (<SAM>) {
		chomp;
		s/^\s+|\s+$//g;
		my ($sam, $pop, $plat) = split(/\t/, $_);
		$sample_pop{$sam} = $pop;	
		$total_sample_cnt_hash{$sam} = 1;
	} 
	my $total_sample_count = keys %total_sample_cnt_hash;
	return $total_sample_count;
}

sub populate_snp_impact_hash {
	my ($y, $identifier) = @_;
	my $snp_func_ref = predict_snp_func($y->{CHROM}, $y->{POS}, $y->{REF}, $y->{ALT}->[0]); ### FIXME, do we worry about non-biallelic snps???
	my %snp_func = %$snp_func_ref;
	my $flag = 0;	
	foreach my $transcript (keys %snp_func) {
		foreach my $string (@{$snp_func{$transcript}}) {
			my ($func, $aa) = split(/:/, $string);
			next if (	$func eq "DOWNSTREAM" ||
						$func eq "UPSTREAM" ||
						$func eq "WITHIN_NON_CODING_GENE" ||
						$func eq "SYNONYMOUS_CODING" || 
						$func eq "INTERGENIC" ||
						$func eq "INTRONIC" ||
						$func eq "NMD_TRANSCRIPT" 
			); ##### NMD_TRANSCRIPT is more an annotation of the transcript than that of a SNP
			print "$transcript-$func\t$aa\t" if ($verbose);	
			my $consequence = "$transcript-$func" . "[" . $aa . "]";
			$snp_impact{$identifier}{$consequence} = 1;
			$flag = 1;
		}
	}
	return $flag;
}        

sub predict_snp_func {
	my ($chr, $pos, $allele1, $allele2) = @_;
	
	#### get Ensembl API registry ####
	my $reg = 'Bio::EnsEMBL::Registry';
	$reg->load_registry_from_db(-host => 'ensembldb.ensembl.org',-user => 'anonymous');
	
	my $vfa	= $reg->get_adaptor('human', 'variation', 'variationfeature');
	my $sa	= $reg->get_adaptor('human', 'core', 'slice');

	my $genotype = $allele1 . "/" . $allele2; ##In SNP predictor, the phase information is not kept.
	#print "Allels are $genotype\n";	
	my $identifier = $chr . ":" . $pos;
	# get a slice for the new feature to be attached to
	my $slice = $sa->fetch_by_region('chromosome', $chr);

	# create a new VariationFeature object
	my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
	  -start => $pos,
	  -end => $pos,
	  -slice => $slice,           		# the variation must be attached to a slice
	  -allele_string => $genotype,    	# the first allele should be the reference allele
	  -strand => 1,						# For 1KG SNPs, use 1
	  -map_weight => 1,
	  -adaptor => $vfa,           		# we must attach a variation feature adaptor
	  -variation_name => $identifier,
	);

	# get the consequence types
	my %tc_string;
	foreach my $con(@{$new_vf->get_all_TranscriptVariations}) {
	  foreach my $string(@{$con->consequence_type}) {
	    my $value = $string . ":" .  $con->pep_allele_string;
	  	push @{$tc_string{$con->transcript->stable_id}},  $value;	  	 
	  }
	}
	return \%tc_string; 
}	

sub post_processing {
	my ($outfile1) = @_;
	`sort -k1nr $outfile1 > $outfile.tmp`;
	`cat $header $outfile.tmp > $outfile.sorted`;
	`rm $header`;
	`rm $outfile.tmp`;
	`rm $outfile1`;
	return 1;
}
	
##################################################################################################################################

=pod

=head1 NAME 

	variant_pattern_finder.pl

=head1 SYNOPSIS

	A script allows one to look for patterns of shared variation between individuals in the same vcf file. 
	To be more specific, in any user-specified chromosomal regions, different samples would have different combination of variations. 
	The finder looks for distinct variation combinations within the region, as well as individuals associated with each variation 
	combination pattern. The finder only focuses on variations that change protein coding sequences such as non-synonymous coding SNPs, 
	splice site changes.

=head1	VERSION

	1.0

=head1	REQUIRED ARGUMENTS

	-vcf		Path to a locally or remotely accessible tabix indexed vcf file. The vcf format is a tab format for presenting variation sites and 
			genotypes data and is described at http://vcftools.sourceforge.net/specs.html. This tool takes both vcf4.0 and vcf4.1 format 
			files.
	-region		Chromosomal region in the format of chr:start-end (1:1000000-100500). As the longer the region is, the more distinctive variant 
			patterns may exist, it is best to work with small regions shorter than several kb.  
	-sample		Path to a tab-delimited file containing sample to population mapping. This information helps to organize the output by population.
			A few lines of example is below:
				
				HG00098 GBR     ILLUMINA
				HG00100 GBR     ILLUMINA
				HG00106 GBR     ILLUMINA

=head1	OPTIONAL ARGUMENTS

	-output_dir		Directory where you want the output file be 
	-print_all		By default, only variants that change a coding sequence would be considered. If -print_all is set, all variants will be considered. 
	-expanded_view		By default -expanded_view is OFF and the finder does not distinguish sites of homozygous reference with those with no data, 
				therefore the number of distinctive combinations of variations is minimized; it offers a simplified and clear variation landscape 
				in the region. The expanded view treats homozygous reference sites and no genotype data sites differently; allows one to see the 
				data with more accuracy. 
	-verbose		print progress along the way				
	-help			Print out help menu

=head1	OUTPUT FILE

	The output file is a tab delimited file. It can be found in the output_dir with this naming convention: chr_start-end.csv.sorted. 
	An example is chr17_41256206-41257206.csv.sorted.  It has two lines of headers.

	Header line 1,		chromosome and chromosomal position of the variation, separated by ":", followed by variation rs number and the reference allele in a 
				square parenthesis.  When rs number is not available, a "." is used instead.
				example: 17:41256702_rs12345[G] or 17:5123456_.[A] 
					
	Header line 2,		functional consequences of the SNP on transcript specified. If multiple transcripts are affected by the variant, they are all 
				listed. The amino acid change is annotated in a square parenthesis at the end.
				example: ENST00000461719-NON_SYNONYMOUS_CODING[T/I]
	
	After the first two header lines, the rest are data content. It has a freq column, a genotype section and a sample section.
	 
	Freq column:		it gives the frequency of the given variant genotype combination in the file
	
	Genotype section:	this section contains individual genotypes. Each column corresponding to one variant.  
				By default, "-" represents genotypes that are either homozygous reference or no data. When -expanded_view is set, "./." represents 
				sites with no genotype data and REF|REF for homozygous reference sites. 
 	
	Sample section:		This section contains samples that have the given pattern of variants. The samples are organized by population. 

=head1 EXAMPLE

perl ~/ReseqTrack/scripts/variation_data/variant_pattern_finder.pl \
-vcf ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20101123/interim_phase1_release/ALL.chr14.phase1.projectConsensus.genotypes.vcf.gz \
-sample /nfs/1000g-archive/vol1/ftp/release/20101123/interim_phase1_release/interim_phase1.20101123.ALL.panel \
-region 14:106329408-106329468 \
-verbose

perl ~/ReseqTrack/scripts/variation_data/variant_pattern_finder.pl \
-vcf /nfs/1000g-archive/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz \
-sample /nfs/1000g-archive/vol1/ftp/release/20100804/20100804.ALL.panel \
-region  17:41256206-41256906
