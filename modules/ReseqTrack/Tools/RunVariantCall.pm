package ReseqTrack::Tools::RunVariantCall;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-reference]   :
      string, path to the reference genome
  Arg [-chrom]   :
      string, optional, chromosome number to be called on
  Arg [-region]   :
      string, optional, chromosomal region to call    
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunVariantCall object.
  Returntype: ReseqTrack::Tools::RunVariantCall
  Exceptions: 
  Example   : my $run_varCall = ReseqTrack::Tools::RunVariantCall->new(
                -program 		=> "variant call program name or path",
                -working_dir 	=> '/path/to/dir/',
                -reference 		=> '/path/to/ref.fa',
				-chrom			=> 20,
				-region			=> 1000000-2000000 );
=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $reference, $chrom, $region)
        = rearrange( [ qw( REFERENCE CHROM REGION ) ], @args);

  $self->reference($reference);
  $self->chrom($chrom);
  $self->region($region);

  return $self;
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

=head2 chrom

  Arg [1]   : ReseqTrack::Tools::RunVariantCall
  Arg [2]   : string, chromosome
  Function  : accessor method for chromosome
  Returntype: string
  Exceptions: n/a
  Example   : my $chrom = $self->chrom;

=cut

sub chrom {
	my ($self, $chr) = @_;
    if ($chr) {
        $self->{'chrom'} = $chr;
    }
    return $self->{'chrom'};
}


=head2 region

  Arg [1]   : ReseqTrack::Tools::RunVariantCall
  Arg [2]   : string, chromosomal region in the format of start-end
  Function  : accessor method for chromosome region
  Returntype: string
  Exceptions: n/a
  Example   : my $region = $self->region;

=cut

sub region {
	my ($self, $region) = @_;
    if ($region) {
        $self->{'region'} = $region;
    }
    return $self->{'region'};
}

=head2 run

  Arg [1]   : ReseqTrack::Tools::RunVariantCall
  Function  : each child object should implement a run method
  Returntype: n/a
  Exceptions: throws as this method should be implement in the child class
  Example   : 

=cut

sub run {
  my ($self) = @_;
  throw(  $self
          . " must implement a run method as ReseqTrack::Tools::RunVariantCall "
          . "does not provide one" );
}

=head2 options							

  Arg [1]   : ReseqTrack::Tools::CallBySamtools or CallByUmake or CallByGATK
  Arg [2]   : string, name of key e.g. "mpileup" 
  Arg [3]   : string, optional, options to be used on command line e.g. "-m 100000000"
  Function  : accessor method for command line options
  Returntype: string, command line options
  Exceptions: n/a
  Example   : my $mpileup_options = $self->options{'mpileup'};

=cut

sub options {
    my ($self, $option_name, $option_value) = @_;

    if (! $self->{'options'}) {
        $self->{'options'} = {};
    }

    throw( "option_name not specified")
        if (! $option_name);

    if ($option_value) {
        $self->{'options'}->{$option_name} = $option_value;
    }

    return $self->{'options'}->{$option_name};
}


=head2 output_raw_bcf

  Arg [1]   : ReseqTrack::Tools::RunVariantCall
  Function  : gets raw bcf file from the list of output files
  Returntype: arrayref of filepaths
  Exceptions: 
  Example   : my $bcf = $self->output_raw_bcf;

=cut

sub output_raw_bcf {
    my $self = shift;
	
    my @output_bcfs;
    foreach my $file (@{$self->output_files}) {
        if ( $file =~ /raw/ ) {
            push( @output_bcfs, $file);
        }
    }

    return \@output_bcfs;
}

=head2 output_flt_vcf

  Arg [1]   : ReseqTrack::Tools::RunVariantCall
  Function  : gets filtered vcf files from the list of output files
  Returntype: arrayref of filepaths
  Exceptions: 
  Example   : my $flt_vcfs = $self->output_flt_vcf;

=cut

sub output_flt_vcf {
    my $self = shift;

    my @output_flt_vcf;
    foreach my $file (@{$self->output_files}) {
        if ( $file !~ /raw/ && $file =~ /vcf/i  ) {
            push @output_flt_vcf, $file;
        }
    }

    return \@output_flt_vcf;
}


=head2 derive_output_file_name 

  Arg [1]   : ReseqTrack::Tools::RunVariantCall::CallByGATK or CallBySamtools or CallByUmake object
  Arg [2]	: algorithm name
  Function  : create an output VCF name based on input file information
  Returntype: file path
  Exceptions: 
  Example   : my $output_file = $self->derive_output_file_name("GATK")->[0];
=cut


sub derive_output_file_name {  
	my ( $self, $algorithm ) = @_;		
	my $first_bam = basename($self->input_files->[0]);
	my @tmp = split(/\./, $first_bam);
	my $first_sample = $tmp[0];
	my $sample_cnt = @{$self->input_files} - 1;
	
	my $output_file;
	if ( $algorithm =~ /gatk/i ) {
		if ($self->region) {
			$output_file = $self->working_dir . "/$first_sample" . "_and_" . $sample_cnt . "_others.chr" . $self->chrom . "_" . $self->region . ".gatk.vcf";
		}
		else {
			$output_file = $self->working_dir . "/$first_sample" . "_and_" . $sample_cnt . "_others.gatk.vcf";
		}
	}
	elsif ( $algorithm =~ /samtools/i ) {
		if ($self->region) {
			$output_file = $self->working_dir . "/$first_sample" . "_and_" . $sample_cnt . "_others.chr" . $self->chrom . "_" . $self->region . ".samtools.raw.bcf";
		}
		else {
			$output_file = $self->working_dir . "/$first_sample" . "_and_" . $sample_cnt . "_others.samtools.raw.bcf";
		}
	}
	elsif ($algorithm =~ /umake/i ) {  ## these are not used, UMAKE has its own naming scheme
		if ($self->region) {
			$output_file = $self->working_dir . "/$first_sample" . "_and_" . $sample_cnt . "_others.chr" . $self->chrom . "_" . $self->region . ".umake.vcf";
		}
		else {
			$output_file = $self->working_dir . "/$first_sample" . "_and_" . $sample_cnt . "_others.umake.vcf";
		}
	}	
	return $self->output_files($output_file);
}	


1;


=pod

=head1 NAME

ReseqTrack::Tools::RunVariantCall

=head1 SYNOPSIS

This is a base class for RunVariantCall objects.
It is a sub class of a ReseqTrack::Tools::RunProgram.
It provides methods that are common to many RunVariantCall child classes 
such as Samtools mpileup and GATK Unified Genotyper.
Child classes should wrap specific alignment algorithms.

=cut
