package ReseqTrack::Tools::RunVariantCall::RunVariantAnnotator;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(basename fileparse);
use ReseqTrack::Tools::FileSystemUtils qw( check_file_exists check_executable);
use ReseqTrack::Tools::RunVariantCall;
use ReseqTrack::Tools::GATKTools;
use File::Path;

@ISA = qw(ReseqTrack::Tools::RunVariantCall ReseqTrack::Tools::GATKTools);

sub new {

  my ( $class, @args ) = @_;	
  my $gatk_obj = $class->ReseqTrack::Tools::GATKTools::new(@args);
  my $variant_call_obj = $class->ReseqTrack::Tools::RunVariantCall::new(@args);
  my $self = {(%$gatk_obj, %$variant_call_obj)};
  bless $self, $class;

  my (	$dbsnp,
  		$input_vcf,
  		$anno,
  		$anno_group,
  		$parameters,
	)
    = rearrange( [ qw( 	DBSNP 
        				INPUT_VCF
    					ANNO
    					ANNO_GROUP
    					parameters		
    				) ], @args);

  
  ### SET DEFAULT

  $self->gatk_path('/nfs/1000g-work/G1K/work/bin/gatk/dist/') if (! $self->gatk_path); #this calls ->program
    
  $self->dbsnp($dbsnp);
  $self->anno($anno);
  $self->anno_group($anno_group);
  $self->parameters($parameters);
  $self->input_vcf($input_vcf);

  return $self;
}


=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::RunVariantAnnotator
  Function  : uses gatk VaraiantAnnotator to annotate VCF files for additional stats,
              which are required by VQSR
              Output is annotated VCF file found in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
    my ($self) = @_; 
         
    my $input_bams = $self->input_files; 
	my $input_vcf_basename = basename($self->input_vcf);
	
    check_file_exists($self->reference);
    check_file_exists(@$input_bams) if ($self->input_files && @{$self->input_files} != 0 );
    check_file_exists($self->input_vcf);  
    check_executable($self->java_exe);

    my $cmd = $self->java_exe . ' ' . $self->jvm_args . ' -jar ';
    $cmd .= $self->gatk_path . '/' . $self->jar_file . " \\\n";

    $cmd .= "-T VariantAnnotator  \\\n";
    $cmd .= "-R " . $self->reference . "  \\\n";
    
    foreach my $bam ( @$input_bams ) { 
        $cmd .= "-I " . $bam . "  \\\n";
    }            
    
    if ( $self->anno &&  @{$self->anno} != 0 ) {
  	  foreach my $anno ( @{$self->anno} ) {
	        $cmd .= "-A " . $anno . "  \\\n";
	    } 
    }
    elsif ( $self->anno_group && @{$self->anno_group} != 0 ) {
  	  foreach my $anno_g ( @{$self->anno_group} ) {
	        $cmd .= "-G " . $anno_g . "  \\\n";
	    } 
    }   
    
    if ( ! $self->anno && ! $self->anno_group && @{$self->anno} == 0 &&  @{$self->anno_group} == 0) {
        throw("Please provide either a list of annotations to be calculated or a annotation_group to be calculated");
    }              
       
    if ( $self->parameters ) {
        foreach my $p ( keys %{$self->parameters} ) {
            $cmd .= "-" . $p . " " . $self->parameters->{$p} . "  \\\n";
        }
    }  
    
    $cmd .= "--variant " . $self->input_vcf . " \\\n";
    
    $cmd .= "--dbsnp " . $self->dbsnp . " \\\n" if ($self->dbsnp);        
    
    $input_vcf_basename =~ s/.gz//g;
    $input_vcf_basename =~ s/.vcf//g;
        
    my $outfile = $self->working_dir . "/" . $input_vcf_basename . ".annotated.vcf"; 
    $cmd .= "-o $outfile \\\n";
    
    if ($self->chrom && $self->region) {
        $cmd .= "-L " . $self->chrom . ":" . $self->region . " \\\n";
    }
    
    print "Running command...........................................\n$cmd\n";
  	
  	$self->output_files($outfile);
  	
  	$self->execute_command_line($cmd);
  	return $self;
}


=head2 input_vcf
  Arg [1]   : ReseqTrack::Tools::RunVariantCall::RunVariantAnnotator
  Arg [2]   : string, required, input_vcf used for run VariantAnnotator
  Function  : accessor method for input_vcf
  Returntype: string
  Exceptions: n/a
  Example   : my $input_vcf = $self->input_vcf;

=cut

sub input_vcf {
  my ($self, $input_vcf) = @_;
  if ($input_vcf) {
    $self->{'input_vcf'} = $input_vcf;
  }
  return $self->{'input_vcf'};
}


=head2 dbsnp
  Arg [1]   : ReseqTrack::Tools::RunVariantCall::RunVariantAnnotator
  Arg [2]   : string, required, dbsnp used for run VariantAnnotator
  Function  : accessor method for dbsnp file
  Returntype: string
  Exceptions: n/a
  Example   : my $dbsnp = $self->dbsnp("/nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/dbsnp_135.b37.vcf");

=cut

sub dbsnp {
  my ($self, $dbsnp) = @_;
  if ($dbsnp) {
    $self->{'dbsnp'} = $dbsnp;
  }
  return $self->{'dbsnp'};
}

=head2 anno

  Arg [1]   : ReseqTrack::Tools::RunVariantAnnotator
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for annotations to be calculated (use the -ls or --list option of VariantAnnotator to see what annotation exist for calculating)
  Returntype: arrayref of strings
  Example   : $self->anno();

=cut

sub anno {
  my ( $self, $arg ) = @_;

  $self->{'anno'} ||= {};
  if ($arg) {
      foreach my $an (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      	    $an =~ s/^\s+|\s+$//g;
      	    $self->{'anno'}->{$an} = 1;	
    	}
  }
  my @annotations = keys %{$self->{'anno'}};  
  return \@annotations;
}


=head2 anno_group

  Arg [1]   : ReseqTrack::Tools::RunVariantAnnotator
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for annotation group to be calculated (use the -ls or --list option of VariantAnnotator to see what anno_grouptation exist for calculating)
  Returntype: arrayref of strings
  Example   : $self->anno_group();

=cut

sub anno_group {
  my ( $self, $arg ) = @_;

  $self->{'anno_group'} ||= {};
  if ($arg) {
      foreach my $an (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      	    $an =~ s/^\s+|\s+$//g;
      	    $self->{'anno_group'}->{$an} = 1;	
    	}
  }
  my @anno_group = keys %{$self->{'anno_group'}};  
  return \@anno_group;
}

1;


=pod


-Xmx2g is too small memory for the chrom20 run
 java -Xmx4g -jar /nfs/1000g-work/G1K/work/bin/gatk/dist/GenomeAnalysisTK.jar \
   -R /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/human_g1k_v37.fasta  \
   -T VariantAnnotator \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_PUR_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_PUR_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_IBS_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CEU_LS454_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CHS_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_ASW_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_IBS_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CLM_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_MXL_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_YRI_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CHB_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_FIN_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_MXL_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_GBR_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CLM_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_LWK_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_YRI_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_JPT_SOLID_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_GBR_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_LWK_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_JPT_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CHB_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CHS_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_ASW_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_FIN_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_TSI_ILLUMINA_chrom20.bam \
   -I /nfs/1000g-work/G1K/scratch/zheng/transpose_bam_chr20/low_coverage_CEU_ILLUMINA_chrom20.bam \
   -o  /nfs/1000g-work/G1K/work/zheng/snp_calling/samtools/all_chr20.new.q_gt_3.copy.annotated.vcf \
   -G StandardAnnotation \
   --variant  /nfs/1000g-work/G1K/work/zheng/snp_calling/samtools/all_chr20.new.q_gt_3.copy.vcf \
   -L 20:1000000-2000000 \
   --dbsnp /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/dbsnp_135.b37.vcf  \

It ran out of memory for the entire chrom20 run, even with 30g memory. Use -L to do it in chunks


#   -A DepthOfCoverage \
#   -A QualByDepth \
#   -A HaplotypeScore \
#   -A MappingQualityRankSumTest \
#   -A ReadPosRankSumTest \


   ## use the -ls or --list option to see what annotation exist foo calculating
   
   Doesn't like *.gz file as input variant file

