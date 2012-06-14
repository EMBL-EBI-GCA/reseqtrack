package ReseqTrack::Tools::RunVariantCall::CallByGATK;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename fileparse);

use base qw(ReseqTrack::Tools::RunVariantCall);
use File::Path;

=head2 new

  Arg [-dbSNP]   :
      string, path to dbsnp vcf file	
  Arg [-dcov]   :
      integer, command line options to use with "samtools mpileup"
  Arg [-stand_emit_conf]   :
      float, The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be emitted (and filtered if less than the calling threshold). 
  Arg [-stand_call_conf]   :
      float, The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be called. 
  Arg [-glm]	:
  	  string, genotype_likelihood_model, choose one of the followings "SNP", "BOTH", "INDEL"
  Arg [-debug_file]		:
		File to print all of the annotated and detailed debugging output.
  Arg [-genotyping_mode]		:
  		String, Default is DISCOVERY. Should we output confident genotypes (i.e. including ref calls) or just the variants? Choose one of the two "DISCOVERY" or "GENOTYPE_GIVEN_ALLELES"
  Arg [-pcr_error_rate]			:
  		Double with default value 1.0E-4
  Arg [-output_mode]			:
  		Should we output confident genotypes (i.e. including ref calls) or just the variants? Default is EMIT_VARIANTS_ONLY
  Arg [-p_nonref_model]			:
  		Non-reference probability calculation model to employ -- EXACT is the default option, while GRID_SEARCH is also available.
  Arg [-minIndelCnt]			:
  		Minimum number of consensus indels required to trigger genotyping run. int with default value 5
  Arg [-min_mapping_quality_score]	:
  		Minimum read mapping quality required to consider a read for calling.
  Arg [-min_base_quality_score]		:
  		Minimum base quality required to consider a base for calling.
  Arg [-metrics_file]	:
  		File to print any relevant callability metrics output.
  Arg [-max_deletion_fraction]	:
  		Maximum fraction of reads with deletions spanning this locus for it to be callable [to disable, set to < 0 or > 1; default:0.05].
  Arg [-indel_heterozygosity]	:
  		Heterozygosity for indel calling. This argument informs the prior probability of having an indel at a site.
  Arg [-heterozygosity]	:
  		Heterozygosity value used to compute prior likelihoods for any locus. 
  Arg [-group]			:
  		One or more classes/groups of annotations to apply to variant calls. Which groups of annotations to add to the output VCF file. 					:
  Arg [-computeSLOD]			:
  		If provided, we will calculate the SLOD. This argument is not enabled by default because it increases the runtime by an appreciable amount.
  Arg [-alleles]				:
  		The set of alleles at which to genotype when in GENOTYPE_MODE = GENOTYPE_GIVEN_ALLELES. 
  Arg [-super_pop_name]	:
      string, default is "unknownPop", used in output file name and collection name
      
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Please see GATK website for detailed descriptions about the paramenters:
  http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html#--output_mode

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallByGATK object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallByGATK
  Exceptions: 
  Example   : my $varCall_byGATK = ReseqTrack::Tools::RunVariantCall::CallByGATK->new(
                -input_files 			=> ['/path/sam1', '/path/sam2'],
                -program 				=> "/path/to/gatk",
                -working_dir 			=> '/path/to/dir/',
                -dbSNP		 			=> '/path/to/dbSNP VCF File/',
                -reference				=> '/path/to/ref/',
                -dcov			 		=> 50,
                -stand_call_conf 		=> 50.0,
                -stand_emit_conf 		=> 10.0,
                -glm					=> 'BOTH',
                -chrom					=> '1',
                -region					=> '1-1000000',
                -output_name_prefix		=> "PHASE1",
                -super_pop_name			=> "ALL"
                 );

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( 	$dbSNP, 
  		$dcov, 
  		$computeSLOD,
  		$stand_emit_conf, 
  		$stand_call_conf, 
  		$debug_file,
  		$glm,
  		$genotyping_mode, 
  		$pcr_error_rate, 
  		$p_nonref_model, 
  		$minIndelCnt,
  		$output_mode,
  		$min_mapping_quality_score,
  		$min_base_quality_score,
  		$metrics_file,
  		$max_deletion_fraction,
  		$indel_heterozygosity,
  		$heterozygosity,
  		$group,
  		$annotation,
  		$alleles,
  		$super_pop_name)
    = rearrange( [ qw( 	DBSNP 
    					DCOV 
    					COMPUTESLOD 
    					STAND_EMIT_CONF 
    					STAND_CALL_CONF 
    					DEBUG_FILE 
    					GLM 
    					GENOTYPING_MODE 
    					PCR_ERROR_RATE 
    					P_NON_REF_MODEL 
    					MININDELCNT 
    					OUTPUT_MODE 
    					MIN_MAPPING_QUALITY_SCORE 
    					MIN_BASE_QUALITY_SCORE 
    					METRICS_FILE 
    					MAX_DELETION_FRACTION 
    					INDEL_HETEROZYGOSITY 
    					HETEROZYGOSITY 
    					GROUP 
    					ANNOTATION 
    					ALLELES
    					SUPER_POP_NAME) ], @args);
  
  my $gatk_path = $self->program;

  ## Set defaults	 
  $self->program("/nfs/1000g-work/G1K/work/bin/gatk/dist/GenomeAnalysisTK.jar") if (! $self->program);
  $self->reference("/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta") if (! $self->reference);
  
  $self->{'options'}->{'dcov'} = 50 unless ($dcov);
  $self->{'options'}->{'stand_emit_conf'} = 10.0 unless ($stand_emit_conf);
  $self->{'options'}->{'stand_call_conf'} = 50.0 unless ($stand_call_conf);
  $self->{'options'}->{'glm'} = 'SNP' unless ($glm);
  
  $self->dbSNP($dbSNP);		
  $self->options('dcov', $dcov);
  $self->options('stand_emit_conf', $stand_emit_conf);
  $self->options('stand_call_conf', $stand_call_conf);
  $self->options('-debug_file', $debug_file);
  $self->options('glm', $glm);
  $self->options('-genotyping_mode', $genotyping_mode);
  $self->options('-pcr_error_rate', $pcr_error_rate);
  $self->options('-output_mode', $output_mode); 
  $self->options('-p_nonref_model', $p_nonref_model); 
  $self->options('minIndelCnt', $minIndelCnt); 
  $self->options('-min_mapping_quality_score', $min_mapping_quality_score); 
  $self->options('-min_base_quality_score', $min_base_quality_score); 
  $self->options('-metrics_file', $metrics_file); 
  $self->options('-max_deletion_fraction', $max_deletion_fraction); 
  $self->options('-indel_heterozygosity', $indel_heterozygosity);
  $self->options('-heterozygosity', $heterozygosity);
  $self->options('-group', $group);
  $self->options('-computeSLOD', $computeSLOD); 
  $self->options('-alleles', $alleles);
  $self->options('-annotation', $annotation);  
  
  if (! $self->dbSNP) {
  	throw("Please provide dbSNP VCF file path for gatk run in parameters\n");
  }
  print "dbSNP is " . $self->dbSNP . "\n";
  
  if (! $self->job_name) {
      $self->generate_job_name;
  }

  $self->super_pop_name($super_pop_name);
  
  return $self;
}

=head2 run

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByGATK
  Function  : uses samtools to call variants in $self->input_files.
  Output is files are stored in $self->output_files
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run {
    my ($self) = @_;  
	
	my $input_bams = $self->input_files; # $input_bams can be an array ref
	my $dir = $self->working_dir;

    my $cmd = "java -jar " . $self->program . " ";
	$cmd .= "-R " . $self->reference . " ";
	
	foreach my $bam ( @$input_bams ) {  ## FIXME: use a list of bam files?
		$cmd .= "-I " . $bam . " ";
	}		
	
	$cmd .= "-T UnifiedGenotyper ";

	my $region;
	if ($self->region) {
		$region = "chr" . $self->chrom . "_" . $self->region; 
	}
	elsif ($self->chrom) {
		$region = "chr" . $self->chrom;
	}	

	my $outfile = $self->derive_output_file_name->[0];
	
	$cmd .= "-o $outfile "; 
	
	$cmd .= "-B:dbsnp,VCF " . $self->dbSNP . " ";
	  
    if ($self->{'options'}) {
        foreach my $tag ( keys ( %{$self->{'options'}} ) ) {
            $cmd .= "-" . $tag . " " . $self->options($tag) . " ";
        }    
    }    
    
    if ($self->chrom) {
        $cmd .= "-L " . $self->chrom;
    	if ($self->region) {
    	    $cmd .= ":" . $self->region;
    	}
    }	        
	
	print "Running command.................................\n";
	my $exit;
	eval{
		$exit = $self->execute_command_line($cmd);
	};
	if( $exit > 0 ){
		throw("Failed to run command\n$cmd\n". @_  . " exit code $exit");
	} 
    
    return $self;
}

=head2 dbSNP

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByGATK
  Arg [2]   : string, path of dbsnp vcf file
  Function  : accessor method for dbsnp
  Returntype: string
  Exceptions: n/a
  Example   : my $dbSNP = $self->dbSNP;

=cut

sub dbSNP {
  my ($self, $dbSNP) = @_;
  if ($dbSNP) {
    $self->{'dbSNP'} = $dbSNP;
  }
  return $self->{'dbSNP'};
}

=head2 super_pop_name
  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByGATK
  Arg [2]   : string, required, super_pop_name used for calling the variants, can be things like EUR, ALL (for all pop) and unknownPop, it will be 
			used in output file names
  Function  : accessor method for super_pop_name
  Returntype: string
  Exceptions: n/a
  Example   : my $super_pop_name = $self->super_pop_name;

=cut

sub super_pop_name {
  my ($self, $super_pop_name) = @_;
  if ($super_pop_name) {
    $self->{'super_pop_name'} = $super_pop_name;
  }
  return $self->{'super_pop_name'};
}

=head2 derive_output_file_name 

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByGATK object
  Arg [2]	: algorithm name
  Function  : create an output VCF name based on input file information
  Returntype: file path
  Exceptions: 
  Example   : my $output_file = $self->derive_output_file_name->[0];

=cut

sub derive_output_file_name {  
	
	
	
	my ( $self ) = @_;		

	my $sample_cnt = @{$self->input_files};
	my $output_file;
	my $output_dir_by_chr;
	my $out_dir = $self->working_dir;
	$out_dir =~ s/\/$//;
	
	if ( $self->chrom ) {
		$output_dir_by_chr = $out_dir . "/chr" . $self->chrom;
	}
	else{
		$output_dir_by_chr = $out_dir;
	}	
	
	mkpath($output_dir_by_chr) unless (-e $output_dir_by_chr);
	
	if ($self->region) {
		$output_file = $output_dir_by_chr . "/" . $self->output_name_prefix . "_" . $self->super_pop_name . "_of_" . $sample_cnt . "bams.chr" . $self->chrom . "_" . $self->region . ".gatk.vcf";

	}
	else {
		$output_file = $output_dir_by_chr . "/" . $self->output_name_prefix . "_" . $self->super_pop_name . "_of_" . $sample_cnt . "bams.gatk.vcf";
	}

	return $self->output_files($output_file);
	
}	

1;

=pod

=head1 NAME

ReseqTrack::Tools::RunVariantCall::CallByGATK

=head1 SYNOPSIS

This is a class for running unified genotyper in the GATK tool kits to call variants
Please refer to GATK website for details of the algorithm:
http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper
It is a sub class of a ReseqTrack::Tools::RunVariantCall.

When run, the module will generate java command like below:
 
> java -jar /nfs/1000g-work/G1K/work/bin/gatk/dist/GenomeAnalysisTK.jar 
-R /nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta 
-I /nfs/1000g-archive/vol1/ftp/data/NA19240/alignment/NA19240.chrom22.LS454.ssaha2.YRI.high_coverage.20100311.bam 
-T UnifiedGenotyper 
-o NA19240.chrom22.LS454.ssaha2.YRI.high_coverage.20100311.gatk.vcf 
-B:dbsnp,VCF /nfs/1000g-work/G1K/scratch/zheng/reference_genomes/dbsnp132_20101103.vcf.gz 
-dcov 50 
-stand_call_conf 50.0 
-stand_emit_conf 10.0  &


