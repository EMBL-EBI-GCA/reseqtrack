package ReseqTrack::Tools::ConvertBam;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;
use File::Path;
use base qw(ReseqTrack::Tools::RunProgram);

=pod

=head1 NAME

ReseqTrack::Tools::ConvertBam;

=head1 SYNOPSIS

class for converting bam into an alternate format.
will currently produce BigWig (.bw) or BedGraph (.bg) or compressed bam (cram and crai)

You need to have the conversion programs in your path:
genomeCoverageBed
bedGraphToBigWig

Or you need to define the convertion program using the -program option.

=head1 Example 1

my $conv = ReseqTrack::Tools::ConvertBam(
                      -input => '/path/to/bam_file'
                      -chromosomes => '/path/to/chromosomes',
                      -output_format => 'bw'
);
$conv->run;
my ($output_file) = @{$conv->output_files};

=head1 Example 2

my $cram = ReseqTrack::Tools::ConvertBam->new (
	-input_files => $bam,
	-output_format => 'cram',
	-options => $parameter_hash,
	-program => $input{program},
	-reference_fasta => $input{reference_fasta},
	-output_file => $output_cram_name,
);

=cut

sub DEFAULT_OPTIONS {
  return { 'split_reads' => 0, };
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $chromosome_file, $format, $output_file, $reference_fasta ) =
    rearrange(
    [qw( CHROMOSOME_FILE OUTPUT_FORMAT OUTPUT_FILE REFERENCE_FASTA )], @args );

  $self->chromosome_file($chromosome_file) if $chromosome_file;
  $self->output_format($format);
  $self->output_file($output_file)         if $output_file;
  $self->reference_fasta($reference_fasta) if $reference_fasta;

  $self->program('genomeCoverageBed') unless $self->program();

  return $self;
}

sub CONVERSION {
  return {
    'bw'   => \&convert_big_wig,
    'bg'   => \&convert_bed_graph,
    'bed'  => \&convert_bed,
    'bb'   => \&convert_big_bed,
    'cram' => \&convert_cram,

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
    "Chromsome file (2 columns, name and length, tab separated) is required")
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

sub convert_cram {
  my ( $self, $input_file, $output_file, $is_intermediate ) = @_;

  my $heap_size;

  if ( $self->options('heap_size') ) {
    $heap_size = '-Xmx' . $self->options('heap_size');
  }
  else {
    $heap_size = "";
  }

  my $executable = $self->program;
  my @cram_cmd = ( 'java', $heap_size, '-jar', $executable );

  if ( $input_file =~ /bam$/ ) {    ### For creating cram files
    push @cram_cmd, "cram";
    push @cram_cmd, '--input-bam-file', $input_file;
    push @cram_cmd, '--output-cram-file', $output_file if ($output_file);
    push @cram_cmd, '--capture-all-tags' if $self->options('capture-all-tags');
    push @cram_cmd, '--ignore-tags', $self->options('ignore-tags')
      if ( $self->options('ignore-tags') );
    push @cram_cmd, '--preserve-read-names'
      if $self->options('preserve-read-names');
    push @cram_cmd, '--lossy-quality-score-spec',
      $self->options('lossy-quality-score-spec')
      if $self->options('lossy-quality-score-spec');
    push @cram_cmd, '-L', $self->options('L') if $self->options('L');
  }
  elsif ( $input_file =~ /cram$/ ) {    ### For indexing cram files
    push @cram_cmd, "index";
    push @cram_cmd, '--input-file', $input_file;
  }

  push @cram_cmd, '--reference-fasta-file', $self->reference_fasta;

  my $run_cram_cmd = join ' ', @cram_cmd;
  $self->do_conversion( $run_cram_cmd, $output_file, $is_intermediate );
}

sub convert_big_wig {
  my ( $self, $input_file, $output_file, $is_intermediate ) = @_;

  my $temp_file = $self->get_temp_dir() . '/temp.bg';

  $self->convert_bed_graph( $input_file, $temp_file, 1 );

  my @big_wig_cmd =
    ( 'bedGraphToBigWig', $temp_file, $self->chromosome_file, $output_file );

  my $bw_cmd = join ' ', @big_wig_cmd;

  $self->do_conversion( $bw_cmd, $output_file, $is_intermediate );
}

sub convert_bed {
  my ( $self, $input_file, $output_file, $is_intermediate ) = @_;

  my @bam_to_bed_cmd = 'bamToBed';

  push @bam_to_bed_cmd, '-split' if ( $self->options('split_reads') );
  push @bam_to_bed_cmd, '-i', $input_file;
  push @bam_to_bed_cmd, '|',  'sort -k1,1 -k2,2n';
  push @bam_to_bed_cmd, '>',  $output_file;

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
  my $output = $self->working_dir;

  if ( $self->output_file ) {
    $output = $self->output_file;
  }
  else {
    $output .= $self->output_name_prefix() if ( $self->output_name_prefix );
    $output .= $self->job_name;
    $output .= '.';
    $output .= $self->output_format;
    print "OUTPUT: $output$/";

  }
  &{ $subs->{$format} }( $self, $input, $output );

}

sub chromosome_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'chromosome_file'} = $arg;
  }
  elsif ( !$self->{'chromosome_file'} ) {
    my $input_file = $self->input_files->[0];
    my $chromosome_file = $self->get_temp_dir() . '/chrom.sizes';
    $self->execute_command_line(
"samtools view -H $input_file | grep \@SQ | cut -f2,3 | sed 's/SN://' | sed 's/LN://' > $chromosome_file"
    );

    $self->{'chromosome_file'} = $chromosome_file;
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

sub output_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'output_file'} = $arg;
  }
  return $self->{'output_file'};
}

sub reference_fasta {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'reference_fasta'} = $arg;
  }
  return $self->{'reference_fasta'};
}

1;
