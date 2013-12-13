package ReseqTrack::Tools::RunVariantCall::CallByGATK;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw( check_file_exists check_executable);
use ReseqTrack::Tools::RunVariantCall;
use ReseqTrack::Tools::GATKTools;
use File::Path;

@ISA = qw(ReseqTrack::Tools::RunVariantCall ReseqTrack::Tools::GATKTools);

=head2 new
  Arg [-parameters]	:
  	  hashref of different parameters to pass on to the program; 
  Arg [-super_pop_name]    :
      string, default is "unknownPop", used in output file name and collection name
      
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Please see GATK website for detailed descriptions about the paramenters:
  http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallByGATK object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallByGATK
  Exceptions: 
  Example   : my $varCall_byGATK = ReseqTrack::Tools::RunVariantCall::CallByGATK->new(
                -input_files			=> ['/path/sam1', '/path/sam2'],
                -program				=> "/path/to/gatk",
                -working_dir			=> '/path/to/dir/',
       			-parameters				=> {'dcov'=>40, 
        									'stand_emit_conf' => 10.0,
     										'stand_call_conf' => 50.0,
     										'glm' => 'SNP',      										
       										'dbSNP'=>'/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/dbsnp132_20101103.vcf.gz'}      
                -reference				=> '/path/to/ref/',
                -chrom					=> '1',
                -region					=> '1-1000000',
                -output_name_prefix		=> "PHASE1",
                -super_pop_name			=> "ALL"
                 );
=cut

sub new {
  my ( $class, @args ) = @_;

  my $gatk_obj = $class->ReseqTrack::Tools::GATKTools::new(@args);
  my $variant_call_obj = $class->ReseqTrack::Tools::RunVariantCall::new(@args);
  my $self = {(%$gatk_obj, %$variant_call_obj)};
  bless $self, $class;

  return $self;
}

# This is to set default parameters
sub DEFAULT_OPTIONS { return {
        'dcov' => 50,
        'stand_emit_conf' => 30.0,
        'stand_call_conf' => 30.0,
        'glm' => 'SNP',
        };
}

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByGATK
  Function  : uses samtools to call variants in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
    my ($self) = @_;  
    my $input_bams = $self->input_files; # $input_bams can be an array ref
    throw "no input bams" if (!@$input_bams);
    $self->check_jar_file_exists;
    check_file_exists($self->reference);
    $self->check_bai_exists();
    check_executable($self->java_exe);

    if (my $alleles_file = $self->options('alleles')) {
      check_file_exists($alleles_file);
      check_file_exists("$alleles_file.tbi");
    }

    my $output_vcf = $self->working_dir .'/'. $self->job_name . '.vcf.gz';
    $output_vcf =~ s{//}{/};

    my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
    push(@cmd_words, $self->gatk_path . '/' . $self->jar_file);
    push(@cmd_words, '-T', 'UnifiedGenotyper');
    push(@cmd_words, '-R', $self->reference);
    
    foreach my $bam ( @$input_bams ) {
        push(@cmd_words, '-I', $bam);
    }        
    
    if (my $region = $self->chrom) {
      if (defined $self->region_start && defined $self->region_end) {
        $region.= ':' . $self->region_start . '-' . $self->region_end;
      }
      push(@cmd_words, '-L', $region);
    }    

    push(@cmd_words, '-o', $output_vcf);
    
    if ( $self->options ) {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;
      }
        foreach my $tag ( keys ( %{$self->options} ) ) {
            my $value = $self->options->{$tag};
            $tag =~ s/^dbSNP/-dbsnp/i;
            push(@cmd_words, "-$tag", $value);
        }
    }        

    my $cmd = join(' ', @cmd_words);
    $self->output_files($output_vcf);
    $self->execute_command_line ($cmd);
    
    return $self;
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
 
> java -jar /nfs/1000g-work/G1K/work/bin/gatk/dist/GenomeAnalysisTK.jar \
-R /nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta \
-I /nfs/1000g-archive/vol1/ftp/data/NA19240/alignment/NA19240.chrom22.LS454.ssaha2.YRI.high_coverage.20100311.bam \
-T UnifiedGenotyper \
-o NA19240.chrom22.LS454.ssaha2.YRI.high_coverage.20100311.gatk.vcf \
--dbsnp /nfs/1000g-work/G1K/scratch/zheng/reference_genomes/dbsnp132_20101103.vcf.gz \
-dcov 50 \
-stand_call_conf 50.0 \
-stand_emit_conf 10.0  &


#About setting -dcov:

It recommend that -dcov is 10 times of average coverage in the samples;  if setting -dcov as 50 in running UnifiedGenotyper, 
then all positions with a coverage over 50 will be downsampled? And for these positions, only 50 mapped reads will be used for variant calling.

For exomes, a straight DP filter shouldn't be used because the relationship between misalignments and depth isn't clear for capture data.

By default the Unified Genotyper downsamples each sample's coverage to no more than 250x (so there will be at most 250 * number_of_samples reads at 
a site). Unless there is a good reason for wanting to change this value, we suggest using this default value especially for exome processing; allowing 
too much coverage will require a lot more memory to run. When running on projects with many samples at low coverage (e.g. 1000 Genomes with 4x 
coverage per sample) we usually lower this value to about 10 times the average coverage: '-dcov 40'.

## About exome studies

All calls made by unified genotyper are to be filtered by VQSR to remove large number of false positives. It is in this filtering step, on-target exome variants 
are selected (adding    --maxGaussians 6 ???)  Cannot find a argument to take bed files.  Anyway look into the Best Practise section for more details:

  http://www.broadinstitute.org/gsa/wiki/index.php/Best_Practice_Variant_Detection_with_the_GATK_v3#Whole_Exome_experiments



