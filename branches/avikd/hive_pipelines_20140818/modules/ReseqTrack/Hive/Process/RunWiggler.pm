package ReseqTrack::Hive::Process::RunWiggler;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use ReseqTrack::Tools::RunWiggler;


sub param_defaults {
  return {
    program_file          => undef,
    output_format         => undef,
    bedGraphToBigWig_path => undef,
    chrom_sizes_file      => undef,
    dedupe                => undef,
    clobber               => undef,
  };
}

sub run {
  my $self = shift @_;
  $self->param_required('bam');
  
  my $bams                   = $self->param_as_array( 'bam' );
  my $output_format          = $self->param( 'output_format' );
  my $bedGraphToBigWig_path  = $self->param( 'bedGraphToBigWig_path' );
  my $chrom_sizes_file       = $self->param( 'chrom_sizes_file' );
  my $dedupe                 = $self->param( 'dedupe' );
  my $clobber                = $self->param( 'clobber' );
  
  foreach my $bam (@$bams) {
    check_file_exists($bam);
  }
  
  my $wiggler = ReseqTrack::Tools::RunWiggler->new (
              -program                => $self->param('program_file'),
              -input_files            => $bams,
              -working_dir            => $self->output_dir,
              -job_name               => $self->job_name,
              -output_format          => $output_format,
              -bedGraphToBigWig_path  => $bedGraphToBigWig_path,
              -chrom_sizes_file       => $chrom_sizes_file,
              -dedupe                 => $dedupe,
              -clobber                => $clobber,
             );
             
  $self->run_program( $wiggler );
  my $output_files = $wiggler->output_files;    
  
  $self->output_param( 'bedgraph', $output_files) if $output_format eq 'bg';  
  $self->output_param( 'bigwig', $output_files)   if $output_format eq 'bw';     
}

1;