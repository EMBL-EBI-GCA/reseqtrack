package ReseqTrack::Tools::RunVariantCall::CallBySamtools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);
use POSIX qw(ceil);

use base qw(ReseqTrack::Tools::RunVariantCall);

=head2 new

  Arg [-bcftools_path]    :
          string, optional, path to bcftools if it cannot be worked out from samtools path
  Arg [-vcfutils_path]    :
          string, optional, path to vcfutils.pl if it cannot be worked out from samtools path         
  Arg [-parameters]   :
      hashref, command line options to use with "samtools mpileup", "bcftools view" and "vcfutils.pl"
      Here is a list of options: samtools mpileup [-EBug] [-C capQcoef]  [-l list] [-M capMapQ] [-Q minBaseQ] [-q minMapQ] 
      Here is a list of options: bcftools view [-AbFGNQSucgv] [-D seqDict] [-l listLoci] [-i gapSNPratio] [-t mutRate] [-p varThres] [-P prior] [-1 nGroup1] [-d minFrac] [-U nPerm] [-X permThres] [-T trioType] 
  Arg [-super_pop_name]    :
      string, default is "unknownPop", used in output file name and collection name
      
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall::CallBySamtools object.
  Returntype: ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Exceptions: 
  Example   : my $varCall_bySamtools = ReseqTrack::Tools::RunVariantCall::CallBySamtools->new(
                -input_files             => ['/path/sam1', '/path/sam2'],
                -program                 => "/path/to/samtools",
                -working_dir             => '/path/to/dir/',  ## this is working dir and output dir
                -bcftools_path           => '/path/to/bcftools/',
                -vcfutils_path           => '/path/to/vcfutils.pl/',
                -reference               => '/path/to/ref/',
                -parameters              => {	'mpileup'=>'-ug', 
                								'bcfview'=>'-bvcg', 
                								'vcfutils'=>'-D 100' }
                -chrom                   => '1',
                -region                  => '1-1000000',
                -output_name_prefix      => PHASE1,
                -super_pop_name          => EUR (or all)
                );

=cut

sub DEFAULT_OPTIONS { return {
        mpileup => '-EDS -e20 -h100 -L250 -o40 -C50 -m1 -F0.002 -d 250 -P ILLUMINA -ug', # mainly the defaults
        bcfview => '-p 0.5 -vcg', #default options
        vcfutils => '-D 10000000 -d 2 -a 2 -Q 10 -w 3 -W 10 -1 1e-4 -2 1e-100 -3 0 -4 1e-4 -e 1e-4', # vcfutils defaults
        depth_of_coverage => undef,
        };
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $bcftools, $vcfutils, $bgzip)
    = rearrange([ qw(
    BCFTOOLS VCFUTILS BGZIP
    )], @args);    
  
  ## Set defaults
  $self->program('samtools') if (! $self->program);
  $self->bcftools($bcftools|| 'bcftools');
  $self->vcfutils($vcfutils|| 'vcfutils.pl');
  $self->bgzip($bgzip|| 'bgzip');

  return $self;
}

=head
sub DEFAULT_OPTIONS { return {
        'mpileup' => '-ug',
        'bcfview' => '-bvcg',
        'vcfutils' => '-d 2',
        };
}
=cut

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallBySamtools
  Function  : uses samtools to call variants in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: 
  Example   : $self->run();

=cut

sub run_program {
    my ($self) = @_;

    throw("do not have have a reference") if !$self->reference;
    check_file_exists($self->reference);
    check_executable($self->bcftools);
    check_executable($self->vcfutils);
    check_executable($self->bgzip);

    throw("Please provide parameters for running mpileup") if !$self->options('mpileup');
    throw("Please provide parameters for running bcfview") if !$self->options('bcfview');
    throw("Please provide parameters for running vcfutils") if !$self->options('vcfutils');

    my $output_vcf = $self->working_dir .'/'. $self->job_name . '.vcf.gz';
    $output_vcf =~ s{//}{/};

    my @cmd_words;
    push(@cmd_words, $self->program, 'mpileup');
    push(@cmd_words, $self->options->{'mpileup'});

    if (my $coverage = $self->options->{'depth_of_coverage'}) {
      push(@cmd_words, '-d', ceil(5.5*$coverage)); # should exceed the -D flag for vcfutils
    }

    push(@cmd_words, '-f', $self->reference);

    if (my $region = $self->chrom) {
      if (defined $self->region_start && defined $self->region_end) {
        $region.= ':' . $self->region_start . '-' . $self->region_end;
      }
      push(@cmd_words, '-r', $region);
    }    
    push(@cmd_words, @{$self->input_files});

    push(@cmd_words, '|', $self->bcftools, 'view');
    push(@cmd_words, $self->options->{'bcfview'});
    push(@cmd_words, '-', '|');

    push(@cmd_words, $self->vcfutils, 'varFilter');
    push(@cmd_words, $self->options->{'vcfutils'});

    if (my $coverage = $self->options->{'depth_of_coverage'}) {
      push(@cmd_words, '-D', ceil(5*$coverage));
    }

    push(@cmd_words, '|', $self->bgzip, '-c');
    push(@cmd_words, '>', $output_vcf);

    my $cmd = join(' ', @cmd_words);
    $self->output_files($output_vcf);
    $self->execute_command_line ($cmd);

    return $self;

}


sub bcftools{
  my ($self, $bcftools) = @_;
  if ($bcftools) {
    $self->{'bcftools'} = $bcftools;
  }
  return $self->{'bcftools'};
}

sub vcfutils{
  my ($self, $vcfutils) = @_;
  if ($vcfutils) {
    $self->{'vcfutils'} = $vcfutils;
  }
  return $self->{'vcfutils'};
}

sub bgzip{
  my ($self, $bgzip) = @_;
  if ($bgzip) {
    $self->{'bgzip'} = $bgzip;
  }
  return $self->{'bgzip'};
}

1;

=pod

=head1 NAME

ReseqTrack::Tools::RunVariantCall::CallBySamtools

=head1 SYNOPSIS

This is a class for running mpileup and BCFtools to call variants
It is a sub class of a ReseqTrack::Tools::RunVariantCall.

It generates and run command like to call SNPs and indels for multiple diploid individuals:

>samtools mpileup -P ILLUMINA -ugf ref.fa *.bam | bcftools view -bcvg - > var.raw.bcf 
>bcftools view var.raw.bcf | vcfutils.pl varFilter -D 2000 > var.flt.vcf

Petr's parameter when doing the 1kg runs:

samtools mpileup -EDS -C50 -m2 -F0.0005 -d 10000 -P ILLUMINA
bcftools view -p 0.99

and with ploidy 1 for chrY and non-pseudosomal regions of male chrX
(1-60000 and 2699521-154931043). Ploidy can be set as an additional
column to -s option of bcftools view, the accepted values are 1 or 2.

For exome studies, extract on-target sites at the end using a tabix command (see mpileup web site use case example)



