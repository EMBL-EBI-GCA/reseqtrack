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
use Net::FTP;
use Vcf;
use File::Copy;
use Bio::EnsEMBL::Variation::Utils::VEP qw(get_all_consequences);

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
		$cache,
		$cache_dir,
		$outfile,
);

my $output_dir = ".";
my $print_all = 0;
my $expanded_view = 0;
my $verbose = 0;
my $host = 'ensembldb.ensembl.org';
my $user = 'anonymous';


&GetOptions( 
	'host=s'		=>\$host,
	'user=s'		=>\$user,
	'vcf=s'			=>\$file,
	'region=s'		=>\$region,
	'sample_panel_file=s'		=>\$sample_panel,
	'output_dir=s'	=>\$output_dir,
	'output_file=s'	=>\$outfile,
	'expanded_view!'	=>\$expanded_view,
	'print_all!'	=>\$print_all,
	'help!'			=>\$help,
	'verbose!'		=>\$verbose,
	'cache!'		=>\$cache,
	'cache_dir=s'	=>\$cache_dir,
);

if ($help) {
	 exec('perldoc', $0);
}

if ($cache) {
	if ( !$cache_dir) {
		die("When run at -cache mode, -cache_dir needs to be defined to specify the location of Ensembl archive files.
see http://www.ensembl.org/info/docs/variation/vep/vep_script.html#cacheopt for details\n");
	}
	unless ( -e $cache_dir ) {
		die("-cache_dir defined $cache_dir does not exist\n");
	}	 
}

	 
die("please provide vcf file, chromosomal region, and sample panel file") if (!$file || !$region || !$sample_panel);

if (!$outfile) {
	$outfile =  "$output_dir/chr" . $region . ".txt";
	$outfile =~ s/:/_/;
}

open(OUT, ">", $outfile) || die("Cannot open output file $outfile");

my $header = "$output_dir/header";
open(HEADER, ">", $header) || die("Cannot open output header file $header");

my $total_sample_cnt = parse_sample_file($sample_panel);

my $vcf	= Vcf->new(file=>$file, region=>$region, print_header=>1); #print_header=>1 allows print sample name rather than column index
$vcf->parse_header();

while (my $x=$vcf->next_data_hash()) {
        #print join ("\t", keys (%$x)) , "\n"; # to get all possible hash keys:
        #FORMAT  QUAL    ID      CHROM   INFO    FILTER  gtypes  REF     ALT     POS 
        # $x->{REF} is string
        # $x->{ALT} is array
        # use 'perldoc Vcf' to see details of Vcf.pm
        
        my $position = $x->{CHROM}.":".$x->{POS};
        
        my $alt = join ("", @{$x->{ALT}});
        print  "SNP: $position " . $x->{ID} . " REF IS: " . $x->{REF} . " ALT IS: " . $alt . "\n" if ($verbose);

		my $identifier = $position . "_" . $x->{ID};
		
		my $is_functional = populate_snp_impact_hash($x, $identifier, $cache, $cache_dir);
		
		next if ( !$is_functional && !$print_all);
		
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
print "Produced ".$outfile."\n";

print "Job done!  Please find the output file in $outfile\n";

################
##### SUBS #####
################
sub parse_sample_file {
	my ($sample_panel) = @_;
	
	my @sample_panel_lines;
	my %total_sample_cnt_hash;
	if ($sample_panel =~ /ftp:\/\/([\w.]+)(\/\S+)/) {
		my $ftp_host = $1;
		my $path = $2;
		
        my $ftp = Net::FTP->new($ftp_host);
        die("Cannot build ftp object, please provide proper ftp path\n") if (!$ftp);
        $ftp->login or die('Cannot login ' , $ftp->message);

        my $sample_panel_content;
        open my $PANEL, '>', \$sample_panel_content;
        $ftp->get($path, $PANEL) or die ("could not open $sample_panel " , $ftp->message);
        $ftp->quit;

        @sample_panel_lines = split(/\n/, $sample_panel_content);
        close $PANEL;
    }
    else {
        open my $FILE, '<', $sample_panel
            or die("cannot open $sample_panel $!");
        @sample_panel_lines = <$FILE>;
        close $FILE;
    }
    
    foreach ( @sample_panel_lines ) {
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
	my ($y, $identifier, $cache, $cache_dir1) = @_;
	
	my $snp_func_ref;
	my $flag = 0;	
	
	## Make sure indels get the correct ensembl-style coordinates	
	my ($s, $e, $allele_string) = vcf_to_ensembl($y); ### FIXME, when Will has appropriate documentation for VEP.pm, 
													  ### will use parse_line function to create VariationFeature objects from VCF line by line.
													  ### variant_pattern.finder.v2.pl uses the VEP functions 
	
	if ( $cache ) {
		$snp_func_ref = predict_snp_func_in_cache($y->{CHROM}, $s, $e, $allele_string, $cache_dir1);
	}
	else {	
		$snp_func_ref = predict_snp_func($y->{CHROM}, $s, $e, $allele_string,); 
	}
	
	my %snp_func = %$snp_func_ref;
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

### To use the cache mode, the ensembl variation archive file needs to be downloaded first at $cache_dir2
### This uses the Bio::EnsEMBL::Variation::Utils::VEP qw(get_all_consequences)
sub predict_snp_func_in_cache {
	my ($chr, $s, $e, $allele_string, $cache_dir2) = @_;

	my ($tr_cache, $rf_cache, $config, %tc_consequence);
		
	# spoof config
	$config = {
		'dir' => $cache_dir2,
		'chunk_size' => 50000,
		'cache_region_size' => 1000000,
		'terms' => 'display',
		'cache' => 1,
		'compress' => 'zcat',
		'quiet' => 1,
	};
		
	my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new_fast({
		start => $s,
		end => $e,
		chr => $chr,
		allele_string => $allele_string,
		strand => 1,
		map_weight => 1,
	});
	
	foreach my $line_hash(@{get_all_consequences($config, [$new_vf], $tr_cache, $rf_cache)}) { ##### Necessary to FIXME? get_all_consequences take a listref of many variantionFeatures objects
																							   ##### consider to input >1 a time			
		unless ( defined($line_hash->{Amino_acids}) ) {
			$line_hash->{Amino_acids} = '-';
		}	
		
		my $value = $line_hash->{Consequence}  . ":" .  $line_hash->{Amino_acids};
		push @{$tc_consequence{$line_hash->{Feature}}}, $value;
=head			
		print join ("\t", (
			#$new_vf->variation_name,
			$line_hash->{Feature},
			$line_hash->{Consequence},
			$line_hash->{Amino_acids}
		));
		print "\n";
=cut	
	}	
	return \%tc_consequence;
}
		
sub predict_snp_func {
	my ($chr, $s, $e, $allele_string) = @_;
		
	#### get Ensembl API registry ####
	my $reg = 'Bio::EnsEMBL::Registry';
	$reg->load_registry_from_db(-host =>$host,-user =>$user);
	
	my $vfa	= $reg->get_adaptor('human', 'variation', 'variationfeature');
	my $sa	= $reg->get_adaptor('human', 'core', 'slice');

	# get a slice for the new feature to be attached to
	my $slice = $sa->fetch_by_region('chromosome', $chr);

	# create a new VariationFeature object
	my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
	  -start => $s,
	  -end => $e,
	  -slice => $slice,           		# the variation must be attached to a slice
	  -allele_string => $allele_string,    	# the first allele should be the reference allele
	  -strand => 1,						# For 1KG SNPs, use 1
	  -map_weight => 1,
	  -adaptor => $vfa,           		# we must attach a variation feature adaptor
#	  -variation_name => $identifier,
	);

	# get the consequence types
	my %tc_string;
	foreach my $con(@{$new_vf->get_all_TranscriptVariations}) {
	  foreach my $string(@{$con->consequence_type}) {
	    unless (defined $con->pep_allele_string ) {
	        $con->pep_allele_string = '-';
	    }    
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
	move($outfile.".sorted", $outfile);
	return 1;
}


#### VCF file records indels differently from ensembl; VCF file only record the start position of the reference allele, 
#### ensembl actually records start and end of an allele
#### This sub is adopted from Will McLaren's parse_vcf function in VEP.pm, it conver the VCF indels coordinates to ensembl convention
#### It also converts non-biallelic alternative alleles to strings separated by '/' rather than ','.  A,G -->A/G
sub vcf_to_ensembl {
	
	my ($z) = @_;
	
	my $chr = $z->{CHROM};
	my $start = $z->{POS};
	my $end = $z->{POS};
	my $ref = $z->{REF}; #string
	my $alt = join('', @{$z->{ALT}});
	 
	# adjust end coord
    $end += (length($ref) - 1);
   
    # find out if any of the alt alleles make this an insertion or a deletion
    my ($is_indel, $is_sub, $ins_count, $total_count);
    foreach my $alt_allele(split /\,/, $alt) {
        $is_indel = 1 if $alt_allele =~ /D|I/;
        $is_indel = 1 if length($alt_allele) != length($ref);
        $is_sub = 1 if length($alt_allele) == length($ref);
        $ins_count++ if length($alt_allele) > length($ref);
        $total_count++;
    }
   
    # multiple alt alleles?
    if($alt =~ /\,/) {
        if($is_indel) {
            my @alts;

            if($alt =~ /D|I/) {
                foreach my $alt_allele(split /\,/, $alt) {
                    # deletion (VCF <4)
                    if($alt_allele =~ /D/) {
                        push @alts, '-';
                    }

                    elsif($alt_allele =~ /I/) {
                        $alt_allele =~ s/^I//g;
                        push @alts, $alt_allele;
                    }
                }
            }

            else {
                $ref = substr($ref, 1) || '-';
                $start++;

                foreach my $alt_allele(split /\,/, $alt) {
                    $alt_allele = substr($alt_allele, 1);
                    $alt_allele = '-' if $alt_allele eq '';
                    push @alts, $alt_allele;
                }
            }

            $alt = join "/", @alts;
        }

        else {
            # for substitutions we just need to replace ',' with '/' in $alt
            $alt =~ s/\,/\//;
        }
    }

    elsif($is_indel) {
        # deletion (VCF <4)
        if($alt =~ /^D/) {
            my $num_deleted = $alt;
            $num_deleted =~ s/\D+//g;
            $end += $num_deleted - 1;
            $alt = "-";
            $ref .= ("N" x ($num_deleted - 1)) unless length($ref) > 1;
        }

        # insertion (VCF <4)
        elsif($alt =~ /^I/) {    
            $ref = '-';
            $alt =~ s/^I//g;
            $start++;
        }

        # insertion or deletion (VCF 4+)
        elsif(substr($ref, 0, 1) eq substr($alt, 0, 1)) { # if the first base is the same for ref and alt, it indicates an indel
            # chop off first base
            $ref = substr($ref, 1) || '-';
            $alt = substr($alt, 1) || '-';

            $start++;
        }
        
        # SV deletion (VCF 4+)
        elsif($alt =~ /DEL/) {
            my $num_deleted = length($ref); #### FIXME, not perfect solution. Wait for Will's code in ensembl release 65
            $end += $num_deleted - 1;
            $alt = "-";
            $ref .= ("N" x ($num_deleted - 1)) unless length($ref) > 1;
        }
        
        elsif($alt =~ /INS/ ) { #### FIXME
        
        }    
    }
	
	my $allele = $ref.'/'.$alt;
	return ($start, $end, $allele);
    
}    
	
##################################################################################################################################
=pod

=head1 NAME 

	variant_pattern_finder.pl

=head1 SYNOPSIS

	A script allows one to look for patterns of shared variation between individuals in the same vcf file. 
	To be more specific, in any user-specified chromosomal regions, different samples would have different combination of variations. 
	The finder looks for distinct variation combinations within the region, as well as individuals associated with each variation 
	combination pattern. The finder can handle SNPs, short indels and structural variations (SV). It uses Ensembl annotations to assess 
	functional consequences of the variants, the assessment method is more mature for SNPs and indels but less so for SVs.  Users
	have the option to output all variants or functional significant variants (non-synonymous coding, frame-shift etc.).    
	
	If the entire input VCF file is phased, the phasing information of any found pattern is accurate. When the input VCF has both phased and
	unphased variants, it is important to output ALL variants in the region using the -print_all option, in order to interpret accurately
	the phase in the found pattern ("/" and "|").  Please see website http://www.broadinstitute.org/gsa/wiki/index.php/Read-backed_phasing_algorithm
	

=head1	VERSION

	1.0

=head1	REQUIRED ARGUMENTS

	-vcf		Path to a locally or remotely accessible tabix indexed vcf file. The vcf format is a tab delimited format for presenting variation sites and 
			genotypes data and is described at http://vcftools.sourceforge.net/specs.html. This tool takes both vcf4.0 and vcf4.1 format 
			files.
	-region		Chromosomal region in the format of chr:start-end (1:1000000-100500). As the longer the region is, the more distinctive variant 
			patterns may exist, it is best to work with small regions shorter than several kb.  
	-sample_panel_file	Path to a tab-delimited file containing sample to population mapping; the file can be either local or remotely accessible. This information 
				helps to organize the output by population.  A few lines of example is below:
				HG00098 GBR     ILLUMINA
				HG00100 GBR     ILLUMINA
				HG00106 GBR     ILLUMINA

=head1	OPTIONAL ARGUMENTS

	-host			Ensembl database host to which you want to connection, default is 'ensembldb.ensembl.org'
	-user			user name to connect to Ensembl database; default is 'anonymous'
	-output_dir		Directory where you want the output file be 
	-output_file		Output file name, if not specified, the script will create output file name based on rule below 
	-print_all		By default, only variants that change a coding sequence would be considered. If -print_all is set, all variants will be considered. 
	-expanded_view		By default -expanded_view is OFF and the finder does not distinguish sites of homozygous reference with those with no data, 
				therefore the number of distinctive combinations of variations is minimized; it offers a simplified and clear variation landscape 
				in the region. The expanded view treats homozygous reference sites and no genotype data sites differently; allows one to see the 
				data with more accuracy. 
	-cache			When database connection is not stable or not desired, the script can be run at -cache mode.  Ensembl variation archive file needs 
				to be downloaded before hand to run in the cache mode, see http://www.ensembl.org/info/docs/variation/vep/vep_script.html#cacheopt for details
	-cache_dir		Directory where the Ensembl archive file is
	-verbose		print progress along the way				
	-help			Print out help menu
				
=head1	OUTPUT FILE

	The output file is a tab delimited file. If user did not specify an output file name, an output file can be found in the output_dir with this naming convention: chr_start-end.txt. 
	An example is chr17_41256206-41257206.txt.  It has two lines of headers.

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
-sample_panel_file /nfs/1000g-archive/vol1/ftp/release/20101123/interim_phase1_release/interim_phase1.20101123.ALL.panel \
-region 14:106329408-106329468 \
-verbose

perl ~/ReseqTrack/scripts/variation_data/variant_pattern_finder.pl \
-vcf /nfs/1000g-archive/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz \
-sample_panel_file /nfs/1000g-archive/vol1/ftp/release/20100804/20100804.ALL.panel \
-region  17:41256206-41256906

perl ~/ReseqTrack/scripts/variation_data/variant_pattern_finder.pl \
-vcf /nfs/1000g-archive/vol1/ftp/release/20110521/ALL.chr14.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.gz \
-sample_panel_file /nfs/1000g-archive/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel \
-region 14:106329408-106329468 \
-verbose  \
-cache \
-cache_dir /nfs/1000g-work/G1K/work/zheng/.vep/homo_sapiens/65 \
-print_all

perl ~/ReseqTrack/scripts/variation_data/variant_pattern_finder.pl \
-vcf /nfs/1000g-archive/vol1/ftp/release/20110521/ALL.chr14.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.gz \
-sample_panel_file ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/phase1_integrated_calls.20101123.ALL.panel \
-region 14:106329408-106329468 \
-verbose \
-cache \
-cache_dir /nfs/1000g-work/G1K/work/zheng/.vep/homo_sapiens/65 \
-print_all