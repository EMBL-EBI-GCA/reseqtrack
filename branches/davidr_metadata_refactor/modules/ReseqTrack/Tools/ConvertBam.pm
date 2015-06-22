package ReseqTrack::Tools::ConvertBam;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use base qw(ReseqTrack::Tools::RunProgram);

=pod

=head1 NAME

ReseqTrack::Tools::ConvertBam;

=head1 SYNOPSIS

class for converting bam into an alternate format.
will currently produce BigWig (.bw) or BedGraph (.bg)

You need to have the conversion programs in your path:
genomeCoverageBed
bedGraphToBigWig

=head1 Example

my $conv = ReseqTrack::Tools::ConvertBam(
                      -input => '/path/to/bam_file'
                      -chromosomes => '/path/to/chromosomes',
                      -format => 'bw'
);
$conv->run;
my ($output_file) = @{$conv->output_files};

=cut

sub DEFAULT_OPTIONS {
  return { 'split_reads' => 0, };
}

sub CONVERSION {
  return {
    'bw' => \&convert_big_wig,
    'bg' => \&convert_bed_graph,
    'bed' => \&convert_bed,
    'bb'  => \&convert_big_bed,

    #to extend: add format and subroutine.
    #each sub routine should take an input and output file name, 
    #and a boolean to indicate whether the output is an intermediate, or the final product 
  };
}

sub get_valid_formats {
  my $subs = CONVERSION();
  return keys %$subs;
}

sub convert_bed_graph {
  my ( $self, $input_file, $output_file, $is_intermediate ) = @_;

  throw(
    "Chromsome file (2 columns, name and length, tab separated) is required" )
    unless ( $self->chromosome_file );

  my @bed_graph_cmd = ( 'genomeCoverageBed', '-bg' );
  my $chromosome_file = $self->chromosome_file;

  push @bed_graph_cmd, '-split' if ( $self->options('split_reads') );
  push @bed_graph_cmd, '-ibam', $input_file;

  push @bed_graph_cmd, '-g', $chromosome_file;

  push @bed_graph_cmd, '>', $output_file;

  my $bg_cmd = join ' ', @bed_graph_cmd;

  $self->do_conversion( $bg_cmd, $output_file, $is_intermediate );
}

sub convert_big_wig {
  my ( $self, $input_file, $output_file, $is_intermediate ) = @_;

  my $temp_file = $self->get_temp_dir() . '/temp.bg';

  $self->convert_bed_graph( $input_file, $temp_file, 1 );

  my $chromosome_file = $self->chromosome_file;

  my @big_wig_cmd =
    ( 'bedGraphToBigWig', $temp_file, $chromosome_file, $output_file );

  my $bw_cmd = join ' ', @big_wig_cmd;

  $self->do_conversion( $bw_cmd, $output_file, $is_intermediate );
}

sub convert_bed {
  my ( $self, $input_file, $output_file, $is_intermediate ) = @_;

  my @bam_to_bed_cmd = 'bamToBed';
  
  push @bam_to_bed_cmd, '-split' if ( $self->options('split_reads') );
  push @bam_to_bed_cmd, '-i', $input_file;
  push @bam_to_bed_cmd, '|',  'sort -k1,1 -k2,2n';
  push @bam_to_bed_cmd, '>', $output_file;

  my $bam_to_bed_cmd = join ' ', @bam_to_bed_cmd;

  $self->do_conversion( $bam_to_bed_cmd, $output_file, $is_intermediate );
}

sub convert_big_bed {
  my ( $self, $input_file, $output_file, $is_intermediate ) = @_;

  my $temp_file = $self->get_temp_dir() . '/temp.bed';

  $self->convert_bed( $input_file, $temp_file, 1 );

  my $chromosome_file = $self->chromosome_file;

  my @big_wig_cmd =
    ( 'bedToBigBed', $temp_file, $chromosome_file, $output_file );

  my $bw_cmd = join ' ', @big_wig_cmd;

  $self->do_conversion( $bw_cmd, $output_file, $is_intermediate );
}

sub do_conversion {
  my ( $self, $cmd, $output_file, $is_intermediate ) = @_;

  if ($is_intermediate) {
    $self->created_files($output_file);
  }
  else {
    $self->output_files($output_file);
  }

  $self->execute_command_line($cmd);
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $chromosome_file, $format ) =
    rearrange( [qw( CHROMOSOME_FILE OUTPUT_FORMAT )], @args );

  $self->chromosome_file($chromosome_file);
  $self->output_format($format);

  $self->program('genomeCoverageBed') unless $self->program();

  return $self;
}

sub run_program {
  my ($self) = @_;

  my $subs   = CONVERSION();
  my $format = $self->output_format();

  throw("Output format is required") unless $format;
  throw("Did not recognise command $format")
    if ( !defined $subs->{$format} );

  my @input_files = @{ $self->input_files };
  throw("ConvertBam works on a single bam ")
    unless ( scalar(@input_files) == 1 );

  my ($input) = ( @{ $self->input_files } );
  my $output = $self->working_dir . '/';
  $output .= $self->output_name_prefix() if ( $self->output_name_prefix );
  $output .= $self->job_name;
  $output .= '.';
  $output .= $self->output_format;

  &{ $subs->{$format} }( $self, $input, $output );

}

sub chromosome_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'chromosome_file'} = $arg;
  }

  return $self->{'chromosome_file'};
}

sub output_format {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'output_format'} = $arg;
  }
  return $self->{'output_format'};
}

1;
