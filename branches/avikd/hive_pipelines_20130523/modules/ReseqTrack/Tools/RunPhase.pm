package ReseqTrack::Tools::RunPhase;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-reference_config]   :
      string, path to the reference config file
  Arg [-chrom]   :
      string, optional, chromosome number to be called on
  Arg [-region]   :
      string, optional, chromosomal region to call    
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunPhase object.
  Returntype: ReseqTrack::Tools::RunPhase
  Exceptions: 
  Example   : my $run_phase = ReseqTrack::Tools::RunPhase->new(
                -program         => "variant call program name or path",
                -working_dir     => '/path/to/dir/',
                -reference_config        => '/path/to/ref_config',
                -chrom            => 20,
                -region            => 1000000-2000000 );

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $reference_config, $chrom, $region)
        = rearrange( [ qw( REFERENCE_CONFIG CHROM REGION ) ], @args);

  $self->reference_config($reference_config);
  $self->chrom($chrom);
  $self->region($region);

  return $self;
}

sub generate_job_name {
    my $self = shift;
    my $job_name = $self->chrom;
    $job_name.='.'. $self->region if($self->region); ## region is optional
    return $self->job_name($job_name);
}

=head2

  Arg [1]   : ReseqTrack::Tools::RunPhase
  Arg [2]   : string, path of reference config file
  Function  : accessor method for reference config file
  Returntype: string
  Exceptions: n/a
  Example   : my $reference = $self->reference_config;
  
=cut

sub reference_config {
    my ($self, $reference_config) = @_;
    if ($reference_config) {
        $self->{'reference_config'} = $reference_config;
    }
    return $self->{'reference_config'};
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

  Arg [1]   : ReseqTrack::Tools::RunPhase
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

  Arg [1]   : ReseqTrack::Tools::RunPhase
  Function  : each child object should implement a run method
  Returntype: n/a
  Exceptions: throws as this method should be implement in the child class
  Example   : 

=cut

sub run_program {
  my ($self) = @_;
  throw(  $self
          . " must implement a run_program method as ReseqTrack::Tools::RunPhase "
          . "does not provide one" );
}

1;


=pod

=head1 NAME

ReseqTrack::Tools::RunPhase

=head1 SYNOPSIS

This is a base class for RunVariantCall objects.
It is a sub class of a ReseqTrack::Tools::RunProgram.
It provides methods that are common to many RunVariantCall child classes 
such as Samtools mpileup and GATK Unified Genotyper.
Child classes should wrap specific alignment algorithms.

=cut
