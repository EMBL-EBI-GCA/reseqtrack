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
                -program         => "variant call program name or path",
                -working_dir     => '/path/to/dir/',
                -reference         => '/path/to/ref.fa',
                -chrom            => 20,
                -region            => 1000000-2000000 );

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

sub generate_job_name {
    my $self = shift;
    my $job_name = $self->chrom .'.'. $self->region;
    return $self->job_name($job_name);
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

=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunVariantCall
  Function  : each child object should implement a run method
  Returntype: n/a
  Exceptions: throws as this method should be implement in the child class
  Example   : 

=cut

sub run_program {
  my ($self) = @_;
  throw(  $self
          . " must implement a run_program method as ReseqTrack::Tools::RunVariantCall "
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

=head2 intermediate_output_file

  Arg [1]   : ReseqTrack::Tools::RunVariantCall
  Function  : stores/gets intermediate output file such as raw vcf
  Returntype: filepaths
  Exceptions: 
  Example   : my $bcf = $self->intermediate_output_file;

=cut

sub intermediate_output_file { ### FIXME: should this be an array?
    my ($self, $file) = @_;
    if ($file) {
        $self->{'intermediate_output_file'} = $file;
    }
    return     $self->{'intermediate_output_file'};
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

