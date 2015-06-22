
package ReseqTrack::Hive::Process::RunHotspot;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunPeakCall::Hotspot;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists );

sub param_defaults {
  return {
     program_file   => undef,
  };
}


sub run {
  my $self = shift @_;
  
  $self->param_required( 'bam' );
  
  my $bams             =  $self->param_as_array( 'bam' );
  my $options          = $self->param('options');

  
  foreach my $bam ( @$bams ) {
    check_file_exists( $bam );
  }
  
  
  my $hotspot = ReseqTrack::Tools::RunPeakCall::Hotspot->new(
					-input_files   => $bams,
					-working_dir   => $self->output_dir,
					-program       => $self->param('program_file'),
					-job_name      => $self->job_name,
                                        -options       => $options,
                    );


  $self->run_program( $hotspot ); 
  $self->output_param( 'hotspot_bed', $hotspot->output_hotspot_bed );
  $self->output_param( 'peak_bed', $hotspot->output_peak_bed ); 
  
}


1;
