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

  my ( $reference, $chrom, $region_start, $region_end)
        = rearrange( [ qw( REFERENCE CHROM REGION_START REGION_END ) ], @args);

  $self->reference($reference);
  $self->chrom($chrom);
  $self->region_start($region_start);
  $self->region_end($region_end);

  return $self;
}

sub generate_job_name {
    my $self = shift;
    my $chrom = $self->chrom;
    if ($self->chrom =~ /\s+/) {   ## when multiple chroms are input, seperated by white space
        $chrom =~ s/\s+/_/g;
    }       
    my $job_name = $chrom;
    if ($self->region_start && $self->region_end) {
      $job_name .= "." . $self->region_start . '-' . $self->region_end;
    }
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

sub region_start {
    my ($self, $region_start) = @_;
    if ($region_start) {
        $self->{'region_start'} = $region_start;
    }
    return $self->{'region_start'};
}

sub region_end {
    my ($self, $region_end) = @_;
    if ($region_end) {
        $self->{'region_end'} = $region_end;
    }
    return $self->{'region_end'};
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

