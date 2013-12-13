
=pod

=head1 NAME

ReseqTrack::Tools::QC::PPQT

=head1 SYNOPSIS

This is a class for running PhantomPeakQualTools to assess the 
quality of a chip seq experiment.

https://code.google.com/p/phantompeakqualtools/

This module only concerns itself with the qc functionality of ppqt,
not the IDR or peak calling functions

=head1 Example

my $ppqt = ReseqTrack::Tools::QC::PPQT->new(
  -program => '/path/to/ppqt/runspp.r',
  -rscript_path => '/path/to/Rscript'
  -input_files => ['my_bam_file.bam'],
  -job_name => 'bob'
);

my $metrics = $ppqt->run;
my $metrics_file = $ppqt->output_files->[0];

=cut

package ReseqTrack::Tools::QC::PPQT;

use strict;
use warnings;
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Tools::GeneralUtils
  qw(execute_system_command execute_pipe_system_command);

use ReseqTrack::Tools::FileSystemUtils
  qw(check_file_exists check_file_does_not_exist);

use base qw(ReseqTrack::Tools::RunProgram);

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);

  my ( $rscript_path, $keep_metrics ) =
    rearrange( [qw( RSCRIPT_PATH KEEP_METRICS)], @args );
  
  $self->rscript_path($rscript_path);
  $self->keep_metrics_file($keep_metrics);
  #setting defaults
  if ( !$self->rscript_path ) {
    $self->rscript_path('Rscript');
  }
  if ( !$self->program ) {
    $self->program('run_spp.R');
  }

  return $self;
}

sub run_ppqt {
  my ($self) = @_;

  my @cmd_args;
  my $job_name    = $self->job_name;
  my $working_dir = $self->working_dir();
  my $param_file  = $working_dir . '/' . $job_name . '.ppqt_metrics';

  check_file_exists( $self->input_files->[0] );
  check_file_does_not_exist($param_file);

  push @cmd_args, $self->rscript_path;
  push @cmd_args, $self->program;

  push @cmd_args, '-c=' . $self->input_files->[0];
  push @cmd_args, '-out=' . $param_file;

  if ( $self->keep_metrics_file ) {
    $self->output_files($param_file);
  }
  else {
    $self->created_files($param_file);
  }

  $self->execute_command_line( join( ' ', @cmd_args ) );

  return $param_file;
}

sub parse_metrics {
  my ( $self, $param_file ) = @_;
  open my $fh, '<', $param_file or throw("could not open parameter file: $!");
  my $line = <$fh>;
  chomp $line;

  my %metrics;
  (
    $metrics{filename},    $metrics{numReads},
    $metrics{estFraglen},  $metrics{corr_estFragLen},
    $metrics{phantomPeak}, $metrics{corr_phantomPeak},
    $metrics{argmin_corr}, $metrics{min_corr},
    $metrics{nsc},         $metrics{rsc},
    $metrics{quality_tag}
  ) = split /\t/, $line;

  my %quality_decode = (
    '-2' => 'veryLow',
    '-1' => 'Low',
    '0'  => 'Medium',
    '1'  => 'High',
    '2'  => 'VeryHigh',
  );

  $metrics{quality_label} = $quality_decode{ $metrics{quality_tag} };

  return [\%metrics];
}

sub run_program {
  my ($self) = @_;

  my $params_file = $self->run_ppqt;
  my $metrics     = $self->parse_metrics($params_file);

  return $metrics;
}

sub rscript_path {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'rscript_path'} = $arg;
  }
  return $self->{'rscript_path'};
}

sub keep_metrics_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'keep_metrics_file'} = $arg;
  }
  return $self->{'keep_metrics_file'};
}

1;
