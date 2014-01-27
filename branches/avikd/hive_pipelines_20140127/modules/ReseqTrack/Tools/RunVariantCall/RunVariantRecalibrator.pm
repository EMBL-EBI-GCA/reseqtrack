package ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw( check_file_exists check_executable);

use base qw(ReseqTrack::Tools::GATKTools);

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
    -program                    => "/path/to/gatk",
    -reference                     => "/nfs/1000g-work/G1K/work/bin/gatk_resource_bundle/human_g1k_v37.fasta",
    -input_files                => ['/path/vcf'],
    -resources                     => hash_ref, {resource_name->resource_parameters resource_path}
    -use_annotation                => ["QD", "HaplotypeScore", "MQRankSum", "ReadPosRankSum", "MQ"],
    -parameters_VR                => hash_ref, {-maxGaussians}->6, {mode}->"BOTH"
    -working_dir                => path_to_output_dir
    -save_files_from_deletion    => 1,
);
  
=cut


sub new {

  my ( $class, @args ) = @_;    
  my $self = $class->SUPER::new(@args);

  my ( $reference, $resources, $annotations)
        = rearrange( [ qw( REFERENCE RESOURCES ANNOTATIONS ) ], @args);
 
  $self->reference($reference);
  $self->resources($resources);
  $self->annotations($annotations);

  return $self;
}

sub DEFAULT_OPTIONS { return {
        'mode' => 'SNP',
        };
}

sub default_annotations {
  my ($self) = @_;
  if ($self->options('mode') eq 'SNP') {
    return [qw(QD HaplotypeScore MQRankSum ReadPosRankSum FS InbreedingCoeff DP)] # MQ????
  }
  if ($self->options('mode') eq 'INDEL') {
    return [qw(DP ReadPosRankSum FS InbreedingCoeff MQRankSum)]
  }
  return [];
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

    my $input_vcf_arr = $self->input_files; 
    throw("expecting one vcf file") if @$input_vcf_arr != 1;
    check_file_exists($input_vcf_arr->[0]);
    check_file_exists($input_vcf_arr->[0] . '.tbi');
    
    $self->check_jar_file_exists;
    check_file_exists($self->reference);
    check_executable($self->java_exe);
    check_file_exists($_) foreach values %{$self->resources};

    my $output_recal_file = $self->working_dir .'/'. $self->job_name . '.recal';
    $output_recal_file =~ s{//}{/}g;
    my $output_tranches_file = $self->working_dir .'/'. $self->job_name . '.tranches';
    $output_tranches_file =~ s{//}{/}g;

    my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
    push(@cmd_words, $self->gatk_path . '/' . $self->jar_file);
    push(@cmd_words, '-T', 'VariantRecalibrator');
    push(@cmd_words, '-R', $self->reference);
    push(@cmd_words, '-input', @$input_vcf_arr);

    while (my ($resource_string, $resource_file) = each %{$self->resources}) {
      push(@cmd_words, "-resource:$resource_string", $resource_file);
    }
    foreach my $annotation (@{$self->annotations}) {
      push(@cmd_words, '-an', $annotation);
    }

    while (my ($tag, $value) = each %{$self->options}) {
        push(@cmd_words, "-$tag", $value);
    }
    push(@cmd_words, '--recal_file', $output_recal_file);
    push(@cmd_words, '--tranches_file', $output_tranches_file);

    my $cmd = join(' ', @cmd_words);
    $self->output_files($output_recal_file);
    $self->output_files($output_tranches_file);
    $self->output_files($output_tranches_file.".pdf");
    $self->execute_command_line ($cmd);

    return;

}

    
=head2 use_annotation

  Arg [1]   : ReseqTrack::Tools::RunVariantRecalibrator
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for annotation flags to be considered in VQSR (VariantRecalibrator).
  Returntype: arrayref of strings
  Example   : $self->use_annotation();

=cut

sub annotations {
  my ( $self, $arg ) = @_;

  $self->{'annotations'} ||= {};
  if ($arg) {
      foreach my $an (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
              $an =~ s/^\s+//;
              $an =~ s/\s+$//g;
              $self->{'annotations'}->{$an} = 1;    
        }
  }
  my @annotations_arr = keys %{$self->{'annotations'}};
  if (!@annotations_arr) {
    return $self->default_annotations;
  }
  return \@annotations_arr;
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
      throw("expecting a hash for resources") if ref($arg) ne 'HASH';
      while (my ($resource_string, $resource_file) = each %$arg) {
        $resource_file =~ s{//}{/}g;
        $self->{'resources'}->{$resource_string} = $resource_file;
      }
  }
  return $self->{'resources'};
}

=head2 reference

  Arg [1]   : ReseqTrack::Tools::RunVariantCall
  Arg [2]   : string, path of reference file
  Function  : accessor method for reference file
  Returntype: string
  Exceptions: n/a
  Example   : my $reference = $self->reference;

=cut

sub reference {
    my ($self, $reference) = @_;
    if ($reference) {
        $self->{'reference'} = $reference;
    }
    return $self->{'reference'};
}

sub output_pdf_file {
    my $self = shift;
    return (grep { /\.pdf$/ } @{ $self->output_files })[0];
}
sub output_recal_file {
    my $self = shift;
    return (grep { /\.recal$/ } @{ $self->output_files })[0];
}
sub output_tranches_file {
    my $self = shift;
    return (grep { /\.tranches$/ } @{ $self->output_files })[0];
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
   





