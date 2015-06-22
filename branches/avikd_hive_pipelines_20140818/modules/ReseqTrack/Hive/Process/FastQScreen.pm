package ReseqTrack::Hive::Process::FastQScreen;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::QC::FastQScreen;

sub param_defaults {
  return {
    program_file       => undef,
    configuration_file => undef,
    bowtie_parameters  => undef,
    subset             => undef,
    keep_text          => undef,
    keep_graph         => undef,
  };
}


sub run {
   my $self = shift @_;
   
   $self->param_required( 'fastq' );
   $self->param_required( 'fastqscreen_conf' );
   my $fastqs       = $self->param_as_array( 'fastq' );
   my $conf_file    = $self->param( 'configuration_file' );
   my $bowtie_param = $self->param( 'bowtie_parameters' );
   my $subset       = $self->param( 'subset' );
   my $keep_text    = $self->param( 'keep_text' );
   my $keep_graph   = $self->param( 'keep_graph' );
   
   
   foreach my $fastq ( @$fastqs ) {
    check_file_exists( $fastq );
  }
  
  my $fastQScreen = ReseqTrack::Tools::QC::FastQScreen->new(
                  -program             => $self->param('program_file'),
                  -input_files         => $fastqs,
                  -configuration_file  => $conf_file,
                  -bowtie_parameters   => $bowtie_param,
                  -subset              => $subset,
                  -keep_text           => $keep_text,
                  -keep_graph          => $keep_graph,
               );
               
  $self->run_program( $fastQScreen );
  $self->output_param( 'fastqscreen', $fastQScreen->output_files ); 
}

1;