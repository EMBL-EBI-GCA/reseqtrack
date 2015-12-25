
package ReseqTrack::Hive::Process::RunMacs2;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunPeakCall::Macs2;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists );

sub param_defaults {
  return {
     program_file   => undef,
     samtools_path  => undef,
     control_files  => undef,
     fragment_size  => undef,
     options        => { 'nomodel' => 1, }, ## Macs2 tools module also setting it 
  };
}


sub run {
  my $self = shift @_;
  
  $self->param_required( 'bam' );
  $self->param_required( 'broad' );
  
  my $bams             =  $self->param_as_array( 'bam' );
  my $samtools_path    =  $self->param( 'samtools_path' );
  my $control_files    =  $self->param( 'control_files' );
  my $fragment_size    =  $self->param( 'fragment_size' );
  my $broad            =  $self->param( 'broad' );
  my $options          = $self->param('options');

  $$options{'broad'}=$broad;

 
  my $fragment_size_stat;

  my @fragment_size_arr   = split /,/, $fragment_size; # ppqt can find multiple peaks. the first is the most likely, so we use that one 
  ( $fragment_size_stat ) = grep{ $_ > 0 } @fragment_size_arr;  
  throw("No fragment size for @$bams") if !$fragment_size_stat;  

  foreach my $bam ( @$bams ) {
    check_file_exists( $bam );
  }
  
  
  my $macs = ReseqTrack::Tools::RunPeakCall::Macs2->new(
					-input_files   => $bams,
					-working_dir   => $self->output_dir,
					-options       => $self->param('options'),
					-program       => $self->param('program_file'),
					-job_name      => $self->job_name,
					-fragment_size => $fragment_size_stat,
					-control_files => $control_files,
					-samtools_path => $samtools_path,
                                        -options       => $options,
                    );


  $self->run_program( $macs ); 
  $self->output_param( 'bed', $macs->bed_file ); 
  $self->output_param( 'bed_xls', $macs->output_bed_xls );
  $self->output_param( 'support_bed', $macs->output_support_bed );
}


1;
