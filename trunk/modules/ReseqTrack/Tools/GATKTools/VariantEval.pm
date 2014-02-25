=pod

=head1 NAME

ReseqTrack::Tools::CallByGATK::VariantEval

=head1 SYNOPSIS

This is a class for running VariantEval in the GATK tool kits to generate a summary report for a VCF file.

Please refer to GATK website for details of the algorithm:
	http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_varianteval_VariantEval.html#--comp
	http://www.broadinstitute.org/gatk/guide/article?id=2361
	
Note that the GATK we installed in house is an old version; the GATK report may be slightly different from what is described in the web site.

When run, the module will generate java command like below:

java -jar /nfs/1000g-work/G1K/work/bin/gatk/dist/GenomeAnalysisTK.jar \
-T VariantEval \
-R /nfs/1000g-work/G1K/work/REFERENCE/aligners_reference/bwa/grc37/human_g1k_v37.fa \
-o merged.no_select.Omni_and_HM.eval.gatkreport \
-eval merged.vcf.gz \
-l INFO \
--evalModule GenotypeConcordance \
-D /nfs/1000g-archive/vol1/ftp/technical/reference/phase2_mapping_resources/ALL.wgs.dbsnp.build135.snps.sites.vcf.gz \
-comp /nfs/1000g-archive/vol1/ftp/technical/working/20131122_broad_omni/Omni25_genotypes_2141_samples.b37.v2.vcf.gz \
-comp /nfs/1000g-archive/vol1/ftp/technical/working/20110322_hapmap3_grch37/hapmap3.3.genotypes.b37_fwd.vcf.gz

=cut

package ReseqTrack::Tools::GATKTools::VariantEval;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw( check_file_exists check_executable);
use base qw(ReseqTrack::Tools::GATKTools);

=head2 new
  Arg [-comps]	:
  	  string, input comparison files that contain datasets you want to compare your evaluation file to. This arg can be used multiple times 
  Arg [-dbsnp]    :
      string, dbSNP file, this is used to classify the evaluation calls into "known" and "novel"
      
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::CallByGATK::VariantEval object.
  Returntype: ReseqTrack::Tools::CallByGATK::VariantEval
  Exceptions: 
  Example   : my $varEval = ReseqTrack::Tools::CallByGATK::VariantEval->new(
                -input_files			=> ['/path/sam1', '/path/sam2'],
                -program				=> "/path/to/gatk",
                -working_dir			=> '/path/to/dir/',
                -reference				=> '/path/to/ref/',
                -comps					=> '/path/to/comparison_file',
                -dbsnp					=> '/path/to/dbsnp_file'
                 );
=cut

sub new {
  my ( $class, @args ) = @_;
  
  	my $self = $class->SUPER::new(@args);

    my ( $comps, $dbsnp)
        = rearrange( [ qw( 	COMPS
        					DBSNP )], @args);

	$self->comps($comps);
	$self->dbsnp($dbsnp);
	return $self;
}

sub run_program {
	my $self = shift;

    throw "no input vcf file to evaluate" if (!$self->input_files);

    if ($self->comps) {
    	check_file_exists($_) foreach (@{$self->comps});
    }
    check_file_exists($self->dbsnp) if $self->dbsnp;
    $self->check_jar_file_exists;
    check_file_exists($self->reference);
    check_executable($self->java_exe);

    my $cmd = $self->java_exe . ' ' . $self->jvm_args . ' -jar ';
    $cmd .= $self->gatk_path . '/' . $self->jar_file;
    $cmd .= " -T VariantEval ";
    $cmd .= "-l INFO ";
	$cmd .= "--evalModule GenotypeConcordance ";
    $cmd .= "-R " . $self->reference . " ";
    $cmd .= "-D " . $self->dbsnp if ($self->dbsnp);
    $cmd .= " -eval " . $self->input_files->[0] . " ";
    
    if ( $self->comps ) {
	    foreach my $comp_vcf ( @{$self->comps} ) {  
    	    $cmd .= "-comp " . $comp_vcf . " ";
    	}
    }	           

    my $outfile = $self->input_files->[0] . ".gatkreport";
    
    $cmd .= "-o $outfile "; 
    
    $self->output_files($outfile);

    $self->execute_command_line($cmd);
    
    return $self;
}

1;


=head2 comps
  Arg [1]   : ReseqTrack::Tools::GATKTools::VariantEval
  Arg [2]   : string or arrayref of strings, this specify the benchmark data sets such as HAPMAP calls and/or OMNI, AFFY genotype calls 
  Function  : accessor method for comps
  Returntype: array ref of string
  Exceptions: n/a
  Example   : my $comps = $self->comps;
=cut


sub comps {
  my ( $self, $arg ) = @_;

  $self->{'comps'} ||= {};
  if ($arg) {
    foreach my $file (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      $file =~ s{//}{/}g;
      $self->{'comps'}->{$file} = 1;
    }
  }

  my @files = keys %{$self->{'comps'}};
  return \@files;
}

=head2 dbsnp
  Arg [1]   : ReseqTrack::Tools::GATKTools::VariantEval
  Arg [2]   : string, this is the dbsnp vcf file path
  Function  : accessor method for dbsnp
  Returntype: string
  Exceptions: n/a
  Example   : my $dbsnp = $self->dbsnp;

=cut

sub dbsnp {
  my ($self, $dbsnp) = @_;
  if ($dbsnp) {
    $self->{'dbsnp'} = $dbsnp;
  }
  return $self->{'dbsnp'};
}





