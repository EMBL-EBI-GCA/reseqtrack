package ReseqTrack::Tools::RunVariantCall::CallByUmake;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);
use File::Find qw(find);
use File::Copy qw (move);

use base qw(ReseqTrack::Tools::RunVariantCall);

=head2 new

  Arg [-parameters]	:
  	  hashref of different parameters to pass on to the program; 
  Arg [-dbSNP]   :
      string, path to prefix of dbsnp vcf rod files    
  Arg [-hm3_prefix]   :
      string, path to prefix of hapmap3 files
  Arg [-indel_prefix]   :
      string, path to prefix of indel vcf files
  Arg [-target_bed]    :
        string,     
                                             
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallByUmake object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallByUmake
  Exceptions: 
  Example   : my $varCall_byGATK = ReseqTrack::Tools::RunVariantCall::CallByUmake->new(
                -input_files             => ['/path/sam1', '/path/sam2'],
                -program                 => "/path/to/gatk",
                -working_dir             => '/path/to/dir/',
                -reference               => '/path/to/ref/',
                -dbSNP                   => 'prefix to dbSNP rod file',
                -hm3_prefix              => 'prefix to hapmap3 files',
                -indel_prefix            => 'prefix to indel vcf files',
                -parameters					=>{ 'offset_off_target'=> 50,
                								'FILTER_MAX_SAMPLE_DP'=>20,
                								'FFILTER_MIN_SAMPLE_DP'=>0.5,
                								'FILTER_MQ'=>3,
                								'FILTER_FLAG'=>0x0704 }
                -chrom                   => '1',
                -region                  => '1-1000000' );

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

 
  ## Set defaults
  $self->program('umake') if (!$self->program); 
  return $self;
}

sub DEFAULT_OPTIONS { return {
    'offset_off_target' => 0,
    'FILTER_MAX_SAMPLE_DP' => 20,
    'FILTER_MIN_SAMPLE_DP' => 1,
    'FILTER_MQ' => 3,
    'FILTER_FLAG' => '0x0704',
    'unit_chunk' => 5000000,
    'LD_Nsnps' => 10000,
    'LD_overlap' => 1000,
        };
}

sub print_bam_index_file {
    my ($self) = @_;    
    my $input_bams = $self->input_files;
    
    my %sample_bams;
    foreach my $bam ( @$input_bams ) {
      my $sample_name = fileparse($bam, qr/\..*/);
      push @{$sample_bams{$sample_name}}, $bam;
    }    
    
    my $bam_index = $self->get_temp_dir . '/' . $self->job_name . '.bam.index';
    open (INDEX, ">", $bam_index) || throw("Cannot open BAM index file $bam_index $!");
    
    while (my ($sample_name, $sample_bams) = each %sample_bams) {
      print INDEX join("\t", $sample_name, 'ALL', @$sample_bams), "\n";
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
    
    my $config = $self->get_temp_dir . '/' . $self->job_name . '.conf';
    open (CONFIG, ">", $config) || throw("Cannot open configuration file $config $!");

    print CONFIG $fixed_text1;
    print CONFIG "UMAKE_ROOT = " . $self->program . "\n";
    #print CONFIG "INPUT_ROOT = " . $self->input_root . "\n";
    print CONFIG "OUTPUT_ROOT = " . $self->get_temp_dir . "\n"; 
    print CONFIG "BAM_INDEX = " . $self->bam_index . "\n";
    print CONFIG "CHRS = " . $self->chrom . "\n" if ($self->chrom);
    
    print CONFIG "OUT_DIR = " . $self->get_temp_dir . "\n";
    print CONFIG "OUT_PREFIX = " . $self->job_name . "\n"; ### This is the prefix for the MAKEFile
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
## (To allow parallel running of multiple jobs, the UNIFORM_TARGET_BED tage is used\
## to take into chromosome chunks, rather than EXOME target bed files\
## For EXOME studies, run as whole genome, then filter on-target calls afterwards\
## using tabix\
## Always set OFFSET_OFF_TARGET = 0 as this is no longer used as EXOME options)\
###############################################################################\n";    

    print CONFIG $fixed_text2;

    if (defined $self->chrom && defined $self->region_start && defined $self->region_end) { 
      my $bed_file = $self->get_temp_dir . '/' . $self->job_name . '.bed';
      open (BED, ">", $bed_file) || throw("Cannot open bed file $bed_file $!");
      print BED join("\t", $self->chrom, $self->region_start, $self->region_end), "\n";
      close BED;

      print CONFIG "WRITE_TARGET_LOCI = TRUE\n"; # FOR TARGETED SEQUENCING ONLY -- Write loci file when performing pileup
      print CONFIG "UNIFORM_TARGET_BED = $bed_file\n";
      print CONFIG "OFFSET_OFF_TARGET = " . $self->options('offset_off_target') . "\n" if $self->options('offset_off_target');
      print CONFIG "MULTIPLE_TARGET_MAP = \n"; # Target per individual : Each line contains [SM_ID] [TARGET_BED]
      print CONFIG "TARGET_DIR = target\n"; # Directory to store target information
      print CONFIG "SAMTOOLS_VIEW_TARGET_ONLY = TRUE\n";  # When performing samtools view, exclude off-target regions (may make command line too long)
    }

    my $fixed_text3 = 
"#\
###############################################################################\
## RESOURCE FILES : Download the full resources for full genome calling\
###############################################################################\n";

    print CONFIG $fixed_text3;
    
    print CONFIG "REF = " . $self->reference . "\n";
    throw("cannot run without option indel_prefix") if !$self->options('indel_prefix');
    throw("cannot run without option dbSNP") if !$self->options('dbSNP');
    throw("cannot run without option hm3_prefix") if !$self->options('hm3_prefix');
    print CONFIG "INDEL_PREFIX = " . $self->options('indel_prefix') . "\n";
    print CONFIG "DBSNP_PREFIX = " . $self->options('dbSNP') . "\n";
    print CONFIG "HM3_PREFIX = " . $self->options('hm3_prefix') . "\n";
    
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
###############################################################################\n";
    print CONFIG $fixed_text5;
    
    print CONFIG "SAMTOOLS_VIEW_FILTER = -q " . $self->options('FILTER_MQ') . " -F " . $self->options('FILTER_FLAG') . "\n"; # samtools view filter (-q by MQ, -F by flag)\n";

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
#############################################################################\n";

    print CONFIG $fixed_text6;
    
    print CONFIG "UNIT_CHUNK = " . $self->options('unit_chunk') . "\n";      # Chunk size of SNP calling : 5Mb is default\
    print CONFIG "LD_NSNPS = " . $self->options('LD_Nsnps') . "\n";          # Chunk size of genotype refinement : 10,000 SNPs\
    print CONFIG "LD_OVERLAP = " . $self->options('LD_overlap') . "\n";        # Overlapping # of SNPs between chunks : 1,000 SNPs\

    my $fixed_text7 =
"RUN_INDEX_FORCE = FALSE   # Regenerate BAM index file even if it exists\
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

    print CONFIG $fixed_text7;
    $self->config($config);
}    


=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByUmake
  Function  : uses umake to call variants in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut


sub run_program {
    my ($self) = @_;  
    
    $self->print_bam_index_file;
    
    $self->print_config_file;

    my $cmd = "perl " . $self->program . "/scripts/umake.pl --conf " .  $self->config . " --snpcall";
    $self->execute_command_line($cmd);
    
    my $cmd2 = "make -f " . $self->get_temp_dir . '/' . $self->job_name . ".Makefile -j 3"; ## FIXME: 3 is a magic number; should take user input
    
    $self->execute_command_line($cmd2);

    my @vcfs;
    find( sub {push(@vcfs, $File::Find::name) if /filtered.vcf.gz$/}, $self->get_temp_dir);
    foreach my $vcf (@vcfs) {
      my $filename = fileparse($vcf);
      $filename =~ s/filtered/umake/;
      my $output_vcf = $self->working_dir . '/'. $self->job_name . '.' . $filename;
      $output_vcf =~ s{//}{/};
      $self->output_files($output_vcf);
      move($vcf, $output_vcf) ;
    }
    
    return $self;
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

1;

=pod

=head1 NAME

ReseqTrack::Tools::RunVariantCall::CallByUmake

=head1 SYNOPSIS

This is a class for running umake to call variants
It is a sub class of a ReseqTrack::Tools::RunVariantCall.

For EXOME studies, as the tag UNIFORM_TARGET_BED has been used for taken in chromosomal regions, 
no tag is left to take in exome bed file. Call as whole genome by parallel running based on chromosome
chunks stored in input_string table and then filter for on-target variants using tabix at the end.  Always set
OFFSET_OFF_TARGET = 0 as this tag is no longer used as an exome tag.


=cut



