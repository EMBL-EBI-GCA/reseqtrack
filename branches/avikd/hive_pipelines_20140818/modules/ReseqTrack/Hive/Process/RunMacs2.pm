
package ReseqTrack::Hive::Process::RunMacs2;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunPeakCall::Macs2;

sub param_defaults {
  return {
     program_file   => undef,
     samtools_path  => undef,
     control_files  => undef,
     fragment_size  => undef,
     options        => {},
  };
}


sub run {
  my $self = shift @_;
  
  $self->param_required( 'bed' );
  my $beds             =  $self->param_as_array( 'bed' );
  my $samtools_path    =  $self->param( 'samtools_path' );
  my $control_files    =  $self->param( 'control_files' );
  my $fragment_size    =  $self->param( 'fragment_size' );
  
  
  foreach my $bed ( @$beds ) {
    check_file_exists( $bed );
  }
  
  
  my $macs = ReseqTrack::Tools::RunPeakCall::Macs2(
					-input_files   => $beds,
					-working_dir   => $self->output_dir,
					-options       => $self->param('options'),
					-program       => $self->param('program_file'),
					-job_name      => $self->job_name,
					-fragment_size => $fragment_size,
					-control_files => $control_files,
					-samtools_path => $samtools_path,
                    );


  $self->run_program( $macs ); 
  $self->output_param( 'bed', $macs->bed_file ); 
  
}


1;
