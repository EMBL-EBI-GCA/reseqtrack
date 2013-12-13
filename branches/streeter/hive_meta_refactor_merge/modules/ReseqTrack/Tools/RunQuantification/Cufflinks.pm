package ReseqTrack::Tools::RunQuantification::Cufflinks;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;

@ISA = qw(ReseqTrack::Tools::RunProgram);

=pod

=head1 NAME

ReseqTrack::Tools::RunQuantification::Cufflinks;

=head1 SYNOPSIS

class for running Cufflinks. Child class of ReseqTrack::Tools::RunProgram

=head1 Example

my $cuff = ReseqTrack::Tools::RunQuantification::Cufflinks(
                      -input => '/path/to/file'
                      -reference => '/path/to/reference',
                      -options => {'threads' => 4},
                      -paired_length => 3000,
                      -read_group_fields => {'ID' => 1, 'LB' => 'my_lib'},
                      -first_read => 1000,
                      -last_read => 2000,
                      );
$cuff->run;
my $output_file_list = $cuff->output_files;

=cut

sub DEFAULT_OPTIONS {
  return {
    'threads'                        => 1,
    'verbose'                        => 0,
    'quantification'                 => 0,
    'guided_assembly'                => 0,
    'mask_file'                      => undef,
    'bias_correction_fasta'          => undef,
    'multi_read_correct'             => 0,
    'mean_fragment_length'           => undef,
    'fragment_length_sd'             => undef,
    'upper_quartile_norm'            => 0,
    'total_hits_norm'                => 0,
    'compatible_hits_norm'           => 0,
    'max-mle-iterations'             => undef,
    'max-bundle-frags'               => undef,
    'no-effective-length-correction' => 0,
    'supress-gsnap-mdz'              => 0,
  };
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  #input files, we will allow single file name or list of file names
  my ($transcript_annotation) =
    rearrange( [qw( TRANSCRIPT_ANNOTATION)], @args );

  $self->transcript_annotation($transcript_annotation);

  #setting defaults
  if ( !$self->program ) {
    if ( $ENV{cufflinks} ) {
      $self->program( $ENV{cufflinks} . '/cufflinks' );
    }
    else {
      $self->program('cufflinks');
    }
  }

  return $self;
}

sub run_program {
  my ($self) = @_;

  my @cmd_args;

  my $out_dir = $self->working_dir() . '/' . $self->job_name;
  $out_dir =~ s!//!/!g;

  push @cmd_args, $self->program;

  print STDERR join "\t", $self->transcript_annotation,
    $self->options('quantification'), $self->options('guided_assembly'), $/;

  throw("Must supply transcript_annotation with these options")
    if ( !$self->transcript_annotation
    && ( $self->options('quantification') || $self->options('guided_assembly') )
    );

  # General Options:
  push @cmd_args, '-o', $self->working_dir() . '/' . $self->job_name;
  push @cmd_args, '-p', $self->options('threads')
    if ( $self->options('threads') );
  push @cmd_args, '--GTF', $self->transcript_annotation
    if ( $self->transcript_annotation && $self->options('quantification') );

#push @cmd_args, '--GTF-guide', $self->transcript_annotation if ( $self->transcript_annotation && $self->options('guided_assembly') );
  push @cmd_args, '--mask-file', $self->options('mask_file')
    if ( $self->options('mask_file') );
  push @cmd_args, '--frag-bias-correct', $self->options('bias_correction_fasta')
    if ( $self->options('bias_correction_fasta') );
  push @cmd_args, '--multi-read-correct'
    if ( $self->options('multi_read_correct') );
  push @cmd_args, '--library_type', $self->library_type
    if ( $self->library_type );

  # Advanced Abundance Estimation Options:
  for my $param_opt (
    qw(frag-len-mean frag-len-std-dev max-mle-iterations max-bundle-frags ))
  {
    push @cmd_args, "--$param_opt", $self->options($param_opt)
      if ( defined $self->options($param_opt) );
  }

  for my $binary_opt (
    qw(upper_quartile_norm total_hits_norm compatible_hits_norm no-effective-length-correction)
    )
  {
    push @cmd_args, "--$binary_opt" if ( $self->options($binary_opt) );
  }

  # Advanced Assembly Options:
  for my $param_opt (
    qw(
    label min-isoform-fraction pre-mrna-fraction max-intron-length junc-alpha
    small-anchor-fraction min-frags-per-transfrag overhang-tolerance max-bundle-length
    min-intron-length trim-3-avgcov-thresh trim-3-dropoff-frac max-multiread-fraction
    overlap-radius )
    )
  {
    push @cmd_args, "--$param_opt", $self->options($param_opt)
      if ( defined $self->options($param_opt) );
  }

  #Advanced Reference Annotation Based Transcript (RABT) Assembly Options:
  if ( $self->options('guided_assembly') ) {
    for my $param_opt (qw(3-overhang-tolerance intron-overhang-tolerance)) {
      push @cmd_args, "--$param_opt", $self->options($param_opt)
        if ( defined $self->options($param_opt) );
    }

    push @cmd_args, '--no-faux-reads' if ( $self->options('no-faux-reads') );

  }

  # Advanced Program Behavior Options:
  if ( $self->options('verbose') ) {
    push @cmd_args, '--verbose';
  }
  else {
    push @cmd_args, '--quiet';
  }
  push @cmd_args, '--no-update-check';

  push @cmd_args, @{ $self->input_files };

  my $cuff_cmd = join( ' ', @cmd_args );

  my @output_files = (
    'genes.fpkm_tracking', 'isoforms.fpkm_tracking',
    'transcripts.gtf',     'skipped.gtf',
  );

  for my $output_file (@output_files) {
    $self->output_files( $out_dir . '/' . $output_file );
  }

  $self->execute_command_line($cuff_cmd);
}

sub transcript_annotation {
  my ( $self, $a ) = @_;
  if ($a) {
    $self->{'transcript_annotation'} = $a;
  }
  return $self->{'transcript_annotation'};
}

sub library_type {
  my $self = shift;
  if (@_) {
    $self->{'library_type'} = (shift) ? 1 : 0;
  }
  return $self->{'library_type'};
}

1;
