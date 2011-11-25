package ReseqTrack::Tools::RunVariantCall::CallByUmake;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename fileparse);

use base qw(ReseqTrack::Tools::RunVariantCall);

=head2 new
  Arg [-dbSNP]   :
      string, path to prefix of dbsnp vcf rod files	
  Arg [-hm3_prefix]   :
      string, path to prefix of hapmap3 files
  Arg [-indel_prefix]   :
      string, path to prefix of indel vcf files
  Arg [-FILTER_MAX_SAMPLE_DP]   :
      integer, argument for filtering
  Arg [-FILTER_MIN_SAMPLE_DP]	:
  	  integer, argument for filtering
  Arg [-offset_off_target]	:
  	  integer,option for EXOME sequencing, extend target by given # of bases	
  Arg [-target_bed]	:
  	  string, 	  	  
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallByUmake object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallByUmake
  Exceptions: 
  Example   : my $varCall_byGATK = ReseqTrack::Tools::RunVariantCall::CallByUmake->new(
                -input_files 			=> ['/path/sam1', '/path/sam2'],
                -program 				=> "/path/to/gatk",
                -working_dir 			=> '/path/to/dir/',
                -reference				=> '/path/to/ref/',
                -dbSNP					=> 'prefix to dbSNP rod file',
                -hm3_prefix				=> 'prefix to hapmap3 files',
                -indel_prefix			=> 'prefix to indel vcf files',
                -offset_off_target		=> 50,
                -FILTER_MAX_SAMPLE_DP	=> 20,
                -FFILTER_MIN_SAMPLE_DP	=> 0.5,
                -chrom					=> '1',
                -region					=> '1-1000000',
                -output_to_working_dir 	=> 1 );

=cut

### FIXME: write to Hyun about what is the best way to pass chromosomal regions; and how to avoid the to run the program using make

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( 	$offset_off_target, 
  		$dbSNP, 
  		$target_bed, 
  		$indel_prefix, 
  		$hm3_prefix, 
  		$FILTER_MAX_SAMPLE_DP, 
  		$FILTER_MIN_SAMPLE_DP)
    = rearrange( [ qw( 	OFFSET_OFF_TARGET 
    					DBSNP 
    					TARGET_BED 
    					INDEL_PREFIX 
    					HM3_PREFIX 
    					FILTER_MAX_SAMPLE_DP 
    					FILTER_MIN_SAMPLE_DP) ], @args);

  ## Set defaults	 
	$self->{'options'}->{'offset_off_target'} = 50 unless ($offset_off_target);
	$self->{'options'}->{'FILTER_MAX_SAMPLE_DP'} = 20 unless ($FILTER_MAX_SAMPLE_DP);
	$self->{'options'}->{'FILTER_MIN_SAMPLE_DP'} = 0.5 unless ($FILTER_MIN_SAMPLE_DP);
 
	$self->reference("/nfs/1000g-work/G1K/work/bin/umake-resources/ref/human.g1k.v37.fa") if ( !$self );  
	$self->dbSNP("/nfs/1000g-work/G1K/work/bin/umake-resources/dbSNP/dbsnp_129_b37.rod") if (!$dbSNP);
	$self->hm3_prefix("/nfs/1000g-work/G1K/work/bin/umake-resources/HapMap3/hapmap3_r3_b37_fwd.consensus.qc.poly") if (!$hm3_prefix);
	$self->indel_prefix("/nfs/1000g-work/G1K/work/bin/umake-resources/indels/1kg.pilot_release.merged.indels.sites.hg19") if (!$indel_prefix);
 
	$self->program("/nfs/1000g-work/G1K/work/bin/umake") if (!$self->program); 
  
	$self->dbSNP($dbSNP) if ($dbSNP);
	$self->indel_prefix($indel_prefix) if ($indel_prefix);
	$self->hm3_prefix($hm3_prefix) if ($hm3_prefix);
	  	
	$self->options('offset_off_target', $offset_off_target);
	$self->options('FILTER_MAX_SAMPLE_DP', $FILTER_MAX_SAMPLE_DP);
	$self->options('FILTER_MIN_SAMPLE_DP', $FILTER_MIN_SAMPLE_DP);
	  
	throw("When run umake, please specify a chromosome") if (!$self->chrom ); 
	  
	if (! $self->job_name) {
		$self->generate_job_name;
	}
	
	return $self;
}

sub print_bam_index_file {
	my ($self) = @_;	
	my $input_bams = $self->input_files;
	
	my %sample_bams;
	my @samples;
	foreach my $bam ( @$input_bams ) {
		my @tmp = split(/\./, basename($bam) );
		my $sample = $tmp[0];
		push @samples, $sample;
		push @{$sample_bams{$sample}}, $bam;
	}	
	
	my $first_sample = $samples[0];
	my $dir = $self->working_dir;
	`mkdir $dir` unless (-e $self->working_dir);
	my $bam_index = $self->working_dir  . "/$first_sample" . "_and_others.bam.index";
	open (INDEX, ">", $bam_index) || throw("Cannot open BAM index file $bam_index");
	
	foreach my $s ( keys %sample_bams ) {
		print INDEX "$s\tALL\t";
		print INDEX join ("\t", @{$sample_bams{$s}} ), "\n";
	}	 
	close INDEX;
	$self->bam_index($bam_index);

	return $self;
}	


sub print_config_file {
	my ($self) = @_;

	my $fixed_text1 = 
"##################################################################\
# UMAKE CONFIGURATION FILE\
# This configuration file contains run-time configuration of\
# UMAKE SNP calling pipeline\
###############################################################################\
## KEY ELEMENTS TO CONFIGURE :\
###############################################################################\n";

	my $config = $self->working_dir . "/test.conf";  ## FIXME, use a better name
	open (CONFIG, ">", $config) || throw("Cannot open configuration file $config\n");

	print CONFIG $fixed_text1;
	print CONFIG "UMAKE_ROOT = " . $self->program . "\n";
	#print CONFIG "INPUT_ROOT = " . $self->input_root . "\n";
	print CONFIG "OUTPUT_ROOT = " . $self->working_dir . "\n";
	print CONFIG "BAM_INDEX = " . $self->bam_index . "\n";
	print CONFIG "CHRS = " . $self->chrom . "\n" if ($self->chrom);
	
	my $out_vcf = $self->working_dir . "/vcfs/chr" . $self->chrom . "/chr" . $self->chrom . ".filtered.vcf.gz" if ($self->chrom);  		
	$self->output_files($out_vcf);																					  
	#### umake can take > 1 chr but i won't allow it here as it doesn't make sense for parallel runs.
	#### FIXME: When the output is Beagle, Thunder or other things, the output directory will be different
	
	my $makeFile_prefix = basename($self->bam_index);
	
	print CONFIG "OUT_DIR = " . $self->working_dir . "\n";
	print CONFIG "OUT_PREFIX = " . $makeFile_prefix . "\n"; ### This is the prefix for the MAKEFile
	print CONFIG "PED_INDEX = test.ped\n" if ($self->chrom && $self->chrom =~ /X/i); # SAMPLE PED FILE (required only for chrX calling)

	my $fixed_text2 = 
"#\
###############################################################################\	
## STEPS TO RUN : COMMENT OUT TO EXCLUDE CERTAIN STEPS\
##   --snpcall, --extract, --beagle, --thunder commands automatically set them\
###############################################################################\
#RUN_INDEX = TRUE        # create BAM index file\
#RUN_PILEUP = TRUE       # create GLF file from BAM\
#RUN_GLFMULTIPLES = TRUE # create unfiltered SNP calls\
#RUN_VCFPILEUP = TRUE    # create PVCF files using vcfPileup and run infoCollector\
#RUN_FILTER = TRUE       # filter SNPs using vcfCooker\
#RUN_SPLIT = TRUE        # split SNPs into chunks for genotype refinement\
#RUN_BEAGLE = TRUE  # BEAGLE - MUST SET AFTER FINISHING PREVIOUS STEPS\
#RUN_SUBSET = TRUE  # SUBSET FOR THUNDER - MAY BE SET WITH BEAGLE STEP TOGETHER\
#RUN_THUNDER = TRUE # THUNDER - MUST SET AFTER FINISHING PREVIOUS STEPS\
#\
###############################################################################\
## OPTIONS FOR GLFEXTRACT (GLFMULTIPLES, VCFPILEUP, FILTER MUST BE TURNED OFF)\
###############################################################################\
#RUN_EXTRACT = TRUE  # Instead of discovering SNPs, extract genotype liklihood in the site of VCF_EXTRACT\
#VCF_EXTRACT = # whole-genome (gzipped and tabixed) .vcf.gz file to extract the site information to genotype (such as 1000 Genomes site list)\
#\
###############################################################################\
## OPTIONS FOR EXOME/TARGETED SEQUENCING : COMMENT OUT IF WHOLE GENOME SEQUENCING\
###############################################################################\n";	
	
	print CONFIG $fixed_text2;

	print CONFIG "WRITE_TARGET_LOCI = TRUE\n"; # FOR TARGETED SEQUENCING ONLY -- Write loci file when performing pileup
	print CONFIG "UNIFORM_TARGET_BED = " . $self->print_target_bed($self->chrom, $self->region) . "\n" if ($self->chrom); # FIXME, please check if this is used correctly.  path for target bed file
	print CONFIG "OFFSET_OFF_TARGET = " . $self->options('offset_off_target') . "\n" if ($self->options('offset_off_target') ); # Extend target by given # of bases
	print CONFIG "MULTIPLE_TARGET_MAP = \n"; # Target per individual : Each line contains [SM_ID] [TARGET_BED]
	print CONFIG "TARGET_DIR = target\n"; # Directory to store target information
	print CONFIG "SAMTOOLS_VIEW_TARGET_ONLY = TRUE\n";  # When performing samtools view, exclude off-target regions (may make command line too long)

	my $fixed_text3 = 
"#\
###############################################################################\
## RESOURCE FILES : Download the full resources for full genome calling\
###############################################################################\n";

	print CONFIG $fixed_text3;
	
	print CONFIG "REF = " . $self->reference . "\n";
	print CONFIG "INDEL_PREFIX = " . $self->indel_prefix . "\n" if ($self->indel_prefix); # indel VCF prefix, example: $(INPUT_ROOT)/data/indels/1kg.pilot_release.merged.indels.sites.hg19
	print CONFIG "DBSNP_PREFIX = " . $self->dbSNP . "\n" if ($self->dbSNP);  ## FIXME, here should be a prefix, need to produce a dbSNP rod map file
	print CONFIG "HM3_PREFIX =  " . $self->hm3_prefix . "\n" if ($self->hm3_prefix); # HapMap3 polymorphic site prefix, example $(INPUT_ROOT)/data/HapMap/hapmap3_r3_b37_fwd.consensus.qc.poly
	
	my $fixed_text4 =
"#\
###############################################################################\
## BINARIES\
###############################################################################\n";	

	print CONFIG $fixed_text4;
	
	print CONFIG "SAMTOOLS_FOR_PILEUP = " . $self->program . "/bin/samtools-hybrid\n";  # for samtools pileup
	print CONFIG "SAMTOOLS_FOR_OTHERS = " . $self->program . "/bin/samtools-hybrid\n";  # for samtools view and calmd
	print CONFIG "GLFMERGE = " . $self->program . "/bin/glfMerge\n";  #glfMerge when multiple BAMs exist per indvidual
	print CONFIG "GLFMULTIPLES = " . $self->program . "/bin/glfMultiples --minMapQuality 0 --minDepth 1 --maxDepth 10000000 --uniformTsTv -smartFilter\n"; # glfMultiples and options
	print CONFIG "GLFEXTRACT = " . $self->program . "/bin/glfExtract\n"; # glfExtract for obtaining VCF for known sites
	print CONFIG "VCFPILEUP = " . $self->program . "/bin/vcfPileup\n"; # vcfPileup to generate rich per-site information
	print CONFIG "INFOCOLLECTOR = " . $self->program . "/bin/infoCollector\n"; # create filtering statistics
	print CONFIG "VCFMERGE = perl " . $self->program . "/scripts/bams2vcfMerge.pl\n"; # merge multiple BAMs separated by chunk of genomes
	print CONFIG "VCFCOOKER = " . $self->program . "/bin/vcfCooker\n";  # vcfCooker for filtering
	print CONFIG "VCFSUMMARY = perl " . $self->program . "/scripts/vcfSummary.pl\n"; # Get summary statistics of discovered site
	print CONFIG "VCFSPLIT = perl " . $self->program . "/scripts/vcfSplit.pl\n"; # split VCF into overlapping chunks for genotype refinement
	print CONFIG "VCFPASTE = perl " . $self->program . "/scripts/vcfPaste.pl\n"; # vcfPaste to generate filtered genotype VCF
	print CONFIG "BEAGLE = java -Xmx4g -jar " .  $self->program . "/ext/beagle.20101226.jar seed=993478 gprobs=true niterations=50 lowmem=true\n"; # BEAGLE BINARY : NEED TO COPY BEAGLE TO $(UMAKE_ROOT)/ext DIRECTORY BEFORE RUNNING PIPELINE
	print CONFIG "VCF2BEAGLE = perl " . $self->program . "/scripts/vcf2Beagle.pl --PL\n";  # convert VCF (with PL tag) into beagle input
	print CONFIG "BEAGLE2VCF = perl " . $self->program . "/scripts/beagle2Vcf.pl\n";  # convert beagle output to VCF
	print CONFIG "THUNDER = " . $self->program . "/bin/thunderVCF -r 30 --phase --dosage --compact --inputPhased\n"; # MaCH/Thunder genotype refinement step
	print CONFIG "LIGATEVCF = perl " . $self->program . "/scripts/ligateVcf.pl\n"; # ligate multiple phased VCFs while resolving the phase between VCFs
	print CONFIG "BGZIP = " . $self->program . "/ext/bgzip\n";  # NEED TO COPY BGZIP TO $(UMAKE_ROOT)/ext DIRECTORY BEFORE RUNNING PIPELINE
	print CONFIG "TABIX = " . $self->program . "/ext/tabix\n";  # NEED TO COPY TABIX TO $(UMAKE_ROOT)/ext DIRECTORY BEFORE RUNNING PIPELINE
	
	my $fixed_text5 =
"#\
###############################################################################\
## ARGUMENT FOR FILTERING\
###############################################################################\
SAMTOOLS_VIEW_FILTER = -q 20 -F 0x0704 # samtools view filter (-q by MQ, -F by flag)\n";

	print CONFIG $fixed_text5;

	print CONFIG "FILTER_MAX_SAMPLE_DP = " . $self->options('FILTER_MAX_SAMPLE_DP') . "\n";  # Max Depth per Sample (20x default) -- will generate FILTER_MAX_TOTAL_DP automatically\
	print CONFIG "FILTER_MIN_SAMPLE_DP = " . $self->options('FILTER_MIN_SAMPLE_DP') . "\n";  # Min Depth per Sample (0.5x defaul) -- will generate FILTER_MIN_TOTAL_DP automatically\
	print CONFIG "FILTER_ARGS = --write-vcf --filter --maxDP \$(FILTER_MAX_TOTAL_DP) --minDP \$(FILTER_MIN_TOTAL_DP) --maxAB 70 --maxSTR 20 --minSTR -20 --winIndel 5 --maxSTZ 5 --minSTZ -5 --maxAOI 5 # arguments for filtering (refer to vcfCooker for details)\n";

	my $fixed_text6 =
"#\
#############################################################################\
## RELATIVE DIRECTORY UNDER OUT_DIR\
#############################################################################\
BAM_GLF_DIR = glfs/bams   # BAM level GLF\
SM_GLF_DIR = glfs/samples # sample level GLF (after glfMerge if necessary)\
VCF_DIR = vcfs            # unfiltered and filtered VCF\
PVCF_DIR = pvcfs          # vcfPileup results\
SPLIT_DIR = split         # chunks split to multiple overlappingpieces\
BEAGLE_DIR = beagle       # beagle output\
THUNDER_DIR = thunder     # MaCH/thunder output\
GLF_INDEX = glfIndex.ped  # glfMultiples/glfExtract index file info\
#\
#############################################################################\
## OTHER OPTIONS\
#############################################################################\
UNIT_CHUNK = 5000000      # Chunk size of SNP calling : 5Mb is default\
LD_NSNPS = 10000          # Chunk size of genotype refinement : 10,000 SNPs\
LD_OVERLAP = 1000         # Overlapping # of SNPs between chinks : 1,000 SNPs\
RUN_INDEX_FORCE = FALSE   # Regenerate BAM index file even if it exists\
MERGE_BEFORE_FILTER = FALSE # Merge across the chromosome before filtering\
NOBAQ_SUBSTRINGS = SOLID  # Avoid BAQ if the BAM file contains the substring\
ASSERT_BAM_EXIST = FALSE  # Check if BAM file exists\
#\
#############################################################################\
## CLUSTER SETTING : CURRENTLY COMPATIBLE WITH MOSIX PLATFORM\
#############################################################################\
MOS_PREFIX =     # PREFIX FOR MOSIX COMMAND (BLANK IF UNUSED)\
MOS_NODES =      # COMMA-SEPARATED LIST OF NODES TO SUBMIT JOBS\
REMOTE_PREFIX =  # REMOTE_PREFIX : Set if cluster node see the directory differently (e.g. /net/mymachine/[original-dir])\n";

	print CONFIG $fixed_text6;
	$self->config($config);
}	

sub print_target_bed {
	my ($self, $chr, $region) = @_;
	my $bed_file = $self->working_dir . "/test.bed";
	open (BED, ">", $bed_file) || throw("Cannot open bed file $bed_file\n");
	print BED "chr$chr\t";
	$region =~ s/-/\t/ if ($region);
	print BED $region . "\n" if ($region);
	return $self->target_bed($bed_file);
}	

sub target_bed {
	my ($self, $bed_file) = @_;
	if ($bed_file) {
		$self->{'target_bed'} = $bed_file;
	}	
	return 	$self->{'target_bed'};
}

sub indel_prefix {
	my ($self, $indel_prefix) = @_;
  if ($indel_prefix) {
    $self->{'indel_prefix'} = $indel_prefix;
  }
  return $self->{'indel_prefix'};
}
	
sub hm3_prefix {
	my ($self, $hm3_prefix) = @_;
  if ($hm3_prefix) {
    $self->{'hm3_prefix'} = $hm3_prefix;
  }
  return $self->{'hm3_prefix'};
}
		
	
sub bam_index {
	my ($self, $bam_index_file) = @_;
  	if ($bam_index_file) {
    	$self->{'bam_index'} = $bam_index_file;
  	}
  	return $self->{'bam_index'};	
}

sub config {
	my ($self, $config) = @_;
  	if ($config) {
    	$self->{'config'} = $config;
  	}
  	return $self->{'config'};	
}

=head2 run

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByUmake
  Function  : uses samtools to call variants in $self->input_files.
  Output is files are stored in $self->output_files
  Returntype: 
  Exceptions: 
  Example   : $self->run();
=cut


sub run {
    my ($self) = @_;  
	
	$self->print_bam_index_file;
	
	$self->print_config_file;
	
	my $cmd = "/nfs/1000g-work/G1K/work/bin/local-perl/bin/perl " . $self->program . "/scripts/umake.pl --conf " .  $self->config . " --snpcall";
	
	print "Running command.............................................................................................\n";
	my $exit;
	eval{
		$exit = $self->execute_command_line($cmd);
	};
	if( $exit > 0 ){
		throw("Failed to run command\n$cmd\n". @_  . " exit code $exit");
	} 
    
    my $cmd2 = "make -f " . $self->bam_index . ".Makefile -j 2"; ## FIXME: 2 is a magic number; should take user input
    
    print "Running command...............................................................................................\n";
	eval{
		$exit = $self->execute_command_line($cmd2);
	};
	if( $exit > 0 ){
		throw("Failed to run command\n$cmd2\n". @_  . " exit code $exit");
	} 
    
    $self->files_to_delete($self->config);
    $self->files_to_delete($self->bam_index);
    $self->files_to_delete($self->bam_index . ".Makefile");
    $self->files_to_delete($self->target_bed)  if ($self->target_bed);
    
	#files stored in hash files_to_delete will be deleted at the end of the run, by destructor DESTROY 
    
    return $self;
}	
	

=head2 dbSNP

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByUmake
  Arg [2]   : string, path to dbsnp vcf files
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


1;

=pod

=head1 NAME

ReseqTrack::Tools::RunVariantCall::CallByUmake

=head1 SYNOPSIS

This is a class for running umake to call variants
It is a sub class of a ReseqTrack::Tools::RunVariantCall.


=cut



