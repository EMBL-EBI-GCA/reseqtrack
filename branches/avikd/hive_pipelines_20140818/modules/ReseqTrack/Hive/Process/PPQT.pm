package ReseqTrack::Hive::Process::PPQT;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::QC::PPQT;
use ReseqTrack::Tools::Exception qw(throw);

sub param_defaults {
  return {
    program_file      => undef,
    rscript_path      => undef,
    samtools          => undef,
    keep_metrics_file => undef,
    keep_plot         => undef,
    keep_rdata        => undef,
    ppqt_no_dups      => undef,
   };
}

sub run {
  my $self = shift @_;

  $self->param_required( 'bam' );
  $self->param_required( 'rscript_path' );
  $self->param_required( 'samtools' );
    
    
  my $bams              = $self->param_as_array( 'bam' );
  my $rscript_path      = $self->param( 'rscript_path' );
  my $samtools_path     = $self->param( 'samtools' );
  my $keep_metrics_file = $self->param( 'keep_metrics_file' );
  my $keep_plot         = $self->param( 'keep_plot' );
  my $keep_rdata        = $self->param( 'keep_rdata' );
  my $ppqt_no_dups      = $self->param( 'ppqt_no_dups' );
    
    
  foreach my $bam (@$bams) {
    check_file_exists($bam);
  }

  my $metrics_generator = ReseqTrack::Tools::QC::PPQT->new(
      -program       => $self->param('program_file'),
      -rscript_path  => $rscript_path,
      -samtools_path => $samtools_path,
      -input_files   => $bams,
      -keep_metrics  => $keep_metrics_file, 
      -keep_plot     => $keep_plot, 
      -keep_rdata    => $keep_rdata, 
      -no_dups       => $ppqt_no_dups, 
      -working_dir   => $self->output_dir,
      -job_name      => $self->job_name,
    );

  $self->run_program( $metrics_generator );
  $self->output_param( 'ppqt', $metrics_generator->output_files->[0] ); 
}

1;