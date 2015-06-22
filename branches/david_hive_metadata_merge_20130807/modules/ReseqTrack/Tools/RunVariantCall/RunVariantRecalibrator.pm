package ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator;

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

=head2 new

  Arg [-resources]   : 
  	  A list of sites for which to apply a prior probability of being correct but which aren't used by the algorithm
  	  String, parameters and path to resources files in the format of resource_name="resource_parameters, resource_path"; 
  	  Here is an example: dbSNP="known=true,training=false,truth=false,prior=8.0 /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/dbsnp_135.b37.vcf"
  	  If multiple resources are needed, use the tag multiple times.
  	  
  Arg [-parameters_VR]   :
      string, parameters to pass to the RunVariantRecalibrator object; if multiple parameters, use this tag multiple times
      
  Arg [-use_annotation]   :  
      string, an annotation that the recalibration will be based on; when multiple, use the tag repetitively. Examples are "QD", "HaplotypeScore", "MQRankSum",
      "ReadPosRankSum" and "MQ"
  
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Please see GATK website for detailed descriptions about inputs to the VariantRecalibrator program:
  http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html#

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator object.
  Returntype: ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator
  Exceptions: 
  Example   : my $object_VR = ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator->new(
	-program					=> "/path/to/gatk",
	-reference 					=> "/nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/human_g1k_v37.fasta",
	-input_files				=> ['/path/vcf'],
	-resources 					=> hash_ref, {resource_name->resource_parameters resource_path}
	-use_annotation				=> ["QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ"],
	-parameters_VR				=> hash_ref, {-maxGaussians}->6, {mode}->"BOTH"
	-working_dir				=> path_to_output_dir
	-save_files_from_deletion	=> 1,
);
  
=cut


sub new {

  my ( $class, @args ) = @_;	
  my $gatk_obj = $class->ReseqTrack::Tools::GATKTools::new(@args);
  my $variant_call_obj = $class->ReseqTrack::Tools::RunVariantCall::new(@args);
  my $self = {(%$gatk_obj, %$variant_call_obj)};
  bless $self, $class;

  my (	$resources,
  		$use_annotation,
  		$parameters_VR,
  		$save_files_from_deletion,
)
    = rearrange( [ qw( 	RESOURCES 
    					USE_ANNOTATION
    					parameters_VR
    					SAVE_FILES_FROM_DELETION
    				) ], @args);
 
  ### SET DEFAULT
  $self->gatk_path('/nfs/1000g-work/G1K/work/bin/gatk/dist/') if (! $self->gatk_path); #this calls ->program
    
  $self->resources($resources);
  $self->use_annotation($use_annotation);
  $self->options($parameters_VR);
  $self->save_files_from_deletion($save_files_from_deletion);

  return $self;
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator
  Function  : VariantRecalibrator is the first step in VQSR variant recalibration.
              The output files from this module are used as input to the ApplyRecalibration module
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
    my ($self) = @_; 

    my $input_vcf = $self->input_files; 
    
    check_file_exists($self->reference);
    check_file_exists(@$input_vcf);
    
    check_executable($self->java_exe);

    my $cmd = $self->java_exe . ' ' . $self->jvm_args . ' -jar ';
    $cmd .= $self->gatk_path . '/' . $self->jar_file . " \\\n";

    $cmd .= "-T VariantRecalibrator  \\\n";
    $cmd .= "-R " . $self->reference . "  \\\n";
    
    foreach my $vcf ( @$input_vcf ) { 
        $cmd .= "-input " . $vcf . "  \\\n";
    }        
    
    foreach my $r ( keys %{$self->resources} ) {
        $cmd .= "-resource:" . $r . "," . $self->resources->{$r} . "  \\\n";
    }
    
    foreach my $annotation ( @{$self->use_annotation} ) {
        $cmd .= "-an " . $annotation . "  \\\n";
    } 
    
    $cmd .= "-recalFile " . 	$self->working_dir . "/output.recal \\\n";
    $cmd .= "-tranchesFile " . 	$self->working_dir . "/output.tranches \\\n";
    $cmd .= "-rscriptFile " . 	$self->working_dir . "/output.plots.R \\\n";
  
    if ( $self->options ) {
        foreach my $p ( keys %{$self->options} ) {
            $cmd .= "-" . $p . " " . $self->options->{$p} . "  \\\n";
        }
    }         
    
    #print "Running command...........................................\n$cmd\n";

  	$self->created_files($self->working_dir . "/output.recal");
  	$self->created_files($self->working_dir . "/output.tranches");
  	$self->created_files($self->working_dir . "/output.plots.R");
	  	
  	$self->execute_command_line($cmd);
  	return $self;
}

    
=head2 use_annotation

  Arg [1]   : ReseqTrack::Tools::RunVariantRecalibrator
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for annotation flags to be considered in VQSR (VariantRecalibrator).
  Returntype: arrayref of strings
  Example   : $self->use_annotation();

=cut

sub use_annotation {
  my ( $self, $arg ) = @_;

  $self->{'use_annotation'} ||= {};
  if ($arg) {
      foreach my $an (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      	    $an =~ s/^\s+|\s+$//g;
      	    $self->{'use_annotation'}->{$an} = 1;	
    	}
  }
  my @annotations = keys %{$self->{'use_annotation'}};  
  return \@annotations;
}

=head2 resources

  Arg [1]   : ReseqTrack::Tools::RunVariantRecalibrator
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for resource files.
  Returntype: arrayref of strings
  Example   : $self->resources('dbsnp=path/to/file');

=cut

sub resources {
  my ( $self, $arg ) = @_;

  $self->{'resources'} ||= {};
  if ($arg) {
      if ( ref($arg) eq 'HASH' ) {    
      	foreach my $resource_name ( keys %$arg ) {
      	    my $resource_file = $arg->{$resource_name};
      	    $resource_file =~ s{//}{/}g;
			$self->{'resources'}->{$resource_name} = $resource_file;
    	}
      }
  }

  my %resources = %{$self->{'resources'}};
  
  return \%resources;
}

1;

=pod  

=head1 Example command line   

 java -Xmx4g -jar /nfs/1000g-work/G1K/work/bin/gatk/dist/GenomeAnalysisTK.jar \
   -T VariantRecalibrator \
   -R /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/human_g1k_v37.fasta \
   -input /nfs/1000g-work/G1K/work/zheng/snp_calling/gatk/gatk_all_chr20.vcf.gz \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/hapmap_3.3.b37.sites.vcf \
   -resource:omni,known=false,training=true,truth=false,prior=12.0 /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/1000G_omni2.5.b37.sites.vcf \
   -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 /nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/dbsnp_135.b37.vcf \
   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ \
   -recalFile ./output.recal \
   -tranchesFile ./output.tranches \
   -rscriptFile ./output.plots.R
   
   When use a chrom chunk worth of vcf file, gatk complaints too few number of sites to be calibrated accuratedly, so use the whole chrom vcf
      -input /nfs/1000g-work/G1K/work/zheng/snp_calling/umake/results3_test/chr10.1000000-2500000.filtered_annotated.vcf \  this one does not contain QD, MQRankSum etc.
 It doesn't like the gz vcf files for resource files. input variant file is ok to be gz  
   





