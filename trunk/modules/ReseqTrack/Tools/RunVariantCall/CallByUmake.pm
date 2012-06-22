package ReseqTrack::Tools::RunVariantCall::CallByUmake;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename fileparse);
use ReseqTrack::Tools::FileSystemUtils;
use File::Copy qw (move);

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
  Arg [-FILTER_MIN_SAMPLE_DP]    :
        integer, argument for filtering
  Arg [-offset_off_target]    :
        integer,option for EXOME sequencing, extend target by given # of bases    
  Arg [-target_bed]    :
        string,     
   Arg [-FILTER_MQ]    : 
        integer, filter output calls by samtools view MQ, default is 3 
  Arg [-FILTER_FLAG]    : 
        string, filter output calls by samtools view flag, deafult is 0x0704
  Arg [-unit_chunk]    :
        integer,    size of chromosomal chunk to call variants on; default is 5Mb 
  Arg [-LD_Nsnps]    :
        integer,   Chunk size of genotype refinement, default is 10,000 SNPs
  Arg [-LD_overlap]    :
        integer, overlapping sizes of chunks, default is 1,000 SNPs
                                             
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallByUmake object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallByUmake
  Exceptions: 
  Example   : my $varCall_byGATK = ReseqTrack::Tools::RunVariantCall::CallByUmake->new(
                -input_files             => ['/path/sam1', '/path/sam2'],
                -program                 => "/path/to/gatk",
                -working_dir             => '/path/to/dir/',
                -reference                => '/path/to/ref/',
                -dbSNP                    => 'prefix to dbSNP rod file',
                -hm3_prefix                => 'prefix to hapmap3 files',
                -indel_prefix            => 'prefix to indel vcf files',
                -offset_off_target        => 50,
                -FILTER_MAX_SAMPLE_DP    => 20,
                -FFILTER_MIN_SAMPLE_DP    => 0.5,
                -FILTER_MQ                => 3,
                -FILTER_FLAG            => 0x0704,
                -chrom                    => '1',
                -region                    => '1-1000000' );

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my (     $offset_off_target, 
          $dbSNP, 
          $target_bed, 
          $indel_prefix, 
          $hm3_prefix, 
          $FILTER_MAX_SAMPLE_DP, 
          $FILTER_MIN_SAMPLE_DP,
          $FILTER_MQ,
          $FILTER_FLAG,
          $unit_chunk,
          $LD_Nsnps,
          $LD_overlap)
    = rearrange( [ qw(     OFFSET_OFF_TARGET 
                        DBSNP 
                        TARGET_BED 
                        INDEL_PREFIX 
                        HM3_PREFIX 
                        FILTER_MAX_SAMPLE_DP 
                        FILTER_MIN_SAMPLE_DP
                        FILTER_MQ
                        FILTER_FLAG
                          UNIT_CHUNK
                          LD_NSNPS
                          LD_OVERLAP ) ], @args);

  ## Set defaults     
      
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
    $self->options('FILTER_MQ', $FILTER_MQ);
    $self->options('FILTER_FLAG', $FILTER_FLAG);          
    $self->options('unit_chunk', $unit_chunk) if ($unit_chunk);
    $self->options('LD_Nsnps', $LD_Nsnps) if ($LD_Nsnps);
    $self->options('LD_overlap', $LD_overlap) if ($LD_overlap);
      
    #throw("When run umake, please specify a chromosome and a region") if (!$self->chrom || !$self->region ); 
     ### FIXME, revive above after test 
    
    return $self;
}

sub DEFAULT_OPTIONS { return {
    'offset_off_target' => 0,
    'FILTER_MAX_SAMPLE_DP' => 20,
    'FILTER_MIN_SAMPLE_DP' => 0.5,
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
    my @samples;
    foreach my $bam ( @$input_bams ) {
        my @tmp = split(/\./, basename($bam) );
        my $sample = $tmp[0];
        push @samples, $sample;
        push @{$sample_bams{$sample}}, $bam;
    }    
    
    my $first_sample = $samples[0];
    my $bam_index = $self->umake_output_dir . "/$first_sample" . "_and_others.bam.index";
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
    
    my $config = $self->umake_output_dir . "/tmp.conf"; 
    $self->created_files($config);
    open (CONFIG, ">", $config) || throw("Cannot open configuration file $config\n");

    print CONFIG $fixed_text1;
    print CONFIG "UMAKE_ROOT = " . $self->program . "\n";
    #print CONFIG "INPUT_ROOT = " . $self->input_root . "\n";
    print CONFIG "OUTPUT_ROOT = " . $self->umake_output_dir . "\n"; 
    print CONFIG "BAM_INDEX = " . $self->bam_index . "\n";
    print CONFIG "CHRS = " . $self->chrom . "\n" if ($self->chrom);
    
    my $makeFile_prefix = basename($self->bam_index);
    
    print CONFIG "OUT_DIR = " . $self->umake_output_dir . "\n";
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
## (To allow parallel running of multiple jobs, the UNIFORM_TARGET_BED tage is used\
## to take into chromosome chunks, rather than EXOME target bed files\
## For EXOME studies, run as whole genome, then filter on-target calls afterwards\
## using tabix\
## Always set OFFSET_OFF_TARGET = 0 as this is no longer used as EXOME options)\
###############################################################################\n";    

    print CONFIG $fixed_text2;

    if ($self->region) { 
        print CONFIG "WRITE_TARGET_LOCI = TRUE\n"; # FOR TARGETED SEQUENCING ONLY -- Write loci file when performing pileup
        print CONFIG "UNIFORM_TARGET_BED = " . $self->print_target_bed($self->chrom, $self->region) . "\n" if ($self->chrom); #path for target bed file
        print CONFIG "OFFSET_OFF_TARGET = " . $self->options('offset_off_target') . "\n" if ($self->options('offset_off_target') ); # Extend target by given # of bases
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

sub print_target_bed {
    my ($self, $chr, $region) = @_;
    my $bed_file = $self->umake_output_dir . "/tmp.bed";
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
    return     $self->{'target_bed'};
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

    $self->umake_output_dir($self->get_temp_dir);
    
    $self->print_bam_index_file;
    
    $self->print_config_file;

    my $cmd = "/nfs/1000g-work/G1K/work/bin/local-perl/bin/perl " . $self->program . "/scripts/umake.pl --conf " .  $self->config . " --snpcall";
    $self->execute_command_line($cmd);
    
    print "Running command.............................................................................................\n";
    
    my $cmd2 = "make -f " . $self->bam_index . ".Makefile -j 10"; ## FIXME: 10 is a magic number; should take user input
    
    print "Running command...............................................................................................\n";
    $self->execute_command_line($cmd2);
    
    $self->sieve_umake_output_dir;
    
    return $self;
}    
    
### umake produces an extensive output directory containing many layers of output files.
### This function is to identify the most relevant output file *filtered.vcf.gz; rename it to a unique name for the chunk
### from a generic name such as chr20.filtered.vcf.gz to chr20.20000000-21000000.filtered.vcf.gz
### prepare for this file to be stored in db and other files to be deleted if keep_intermediate_file == 0
sub sieve_umake_output_dir {
    my ($self) = @_;
    my @outfiles_to_del;

    my ($output_files, $outdir_hash) = list_files_in_dir($self->umake_output_dir);

    foreach my $sub_dir ( keys %$outdir_hash ) {
        foreach my $out_file ( @{$outdir_hash->{$sub_dir}} ) {
            my $out_file_path = $sub_dir . "/" . $out_file;
            if ($out_file =~ /filtered.vcf.gz$|filtered.vcf.gz.tbi$/) {
                my $chunk = $self->region;
                my $edited_out_file = $out_file;
                $edited_out_file =~ s/filtered/$chunk.filtered/;  
                $self->output_files($edited_out_file);
                move($out_file_path, $edited_out_file)
                    or throw("Failed to change the output filtered vcf file name $out_file_path: $!");
            }
        }
    }
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

sub umake_output_dir {
    my ($self, $umake_output_dir) = @_;
      if ($umake_output_dir) {
        $self->{'umake_output_dir'} = $umake_output_dir;
      }
      return $self->{'umake_output_dir'};    
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



