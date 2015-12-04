package ReseqTrack::Hive::Process::RunWiggler;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunWiggler;
use ReseqTrack::Tools::FileSystemUtils qw( check_file_exists );

sub param_defaults {
  return {
     program_file            => undef,
     fragment_size           => undef,
     samtools_path           => undef,
     collection_name         => undef,
     bam_type                => undef,
     dedup                   => 0,
     fragment_size_stat_name => 'estFraglen',
     output_format           => 'bg',
     bedGraph_to_bigWig_path => undef,
     chrom_sizes_file        => undef,
     chrom_fasta_file        => undef,
     mappability_tracks      => undef, 
     options                 => {},  
  };
}


sub run {
  my $self = shift @_;
  
  $self->param_required( 'bam' );
  my $bams                    = $self->param_as_array( 'bam' );
  my $samtools_path           = $self->param( 'samtools_path' );
  my $fragment_size           = $self->param( 'fragment_size' );
  my $options                 = $self->param( 'options' );
  my $collection_name         = $self->param( 'collection_name' );
  my $bam_type                = $self->param( 'bam_type' );
  my $fragment_size_stat_name = $self->param( 'fragment_size_stat_name' );
  my $bedGraph_to_bigWig_path = $self->param( 'bedGraph_to_bigWig_path' );
  my $output_format           = $self->param( 'output_format' );
  my $chrom_sizes_file        = $self->param( 'chrom_sizes_file' );
  my $dedup                   = $self->param( 'dedup' );
  my $chrom_fastq_file        = $self->param_required( 'chrom_fasta_file' ); 
  my $mappability_tracks      = $self->param_required( 'mappability_tracks' );
  my $mcr_root                = $self->param_required( 'mcr_root' );
 

  throw( "unknown output format $output_format" ) 
       unless ( $output_format eq 'bg' or $output_format eq 'bw' );

  $fragment_size = 1 unless $fragment_size; ## fix for Dnase seq, need to fix this 

  unless( $fragment_size ){
    throw( 'need collection name and BAM type for fetching fragment size for '.$$bams[0]) 
         if ( !$collection_name or !$bam_type );
  
    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param_required('reseqtrack_db')});
    my $ca = $db->get_CollectionAdaptor;
    my $collection = $ca->fetch_by_name_and_type( $collection_name, $bam_type );

    throw( "failed to fetch any collection for $collection_name and type $bam_type" ) 
         unless $collection;
 
    my $stats = $collection->attributes();
  
    for my $stat (@$stats) {
      if ( $stat->attribute_name eq $fragment_size_stat_name ) { 
          $fragment_size = $stat->attribute_value;
      }
    }  
    $db->dbc->disconnect_when_inactive(1);    
  } 

  my $fragment_size_stat; 
  my @fragment_size_arr = split /,/, $fragment_size; # ppqt can find multiple peaks. the first is the most likely, so we use that one 
  
  foreach my $frag (@fragment_size_arr){
    if( $frag > 0){
      $fragment_size_stat = $frag;
      last;
    }
  }

  throw("No fragment size for @$bams") if !$fragment_size_stat;   
  $fragment_size = $fragment_size_stat;
 
  foreach my $bam ( @$bams ) {
    check_file_exists( $bam );
  }
  
  $$options{'s'}=$chrom_fastq_file;
  $$options{'u'}=$mappability_tracks;
  $$options{'l'}=$fragment_size;
 
  my $wiggler = ReseqTrack::Tools::RunWiggler->new(
			                -input_files           => $bams,
					-working_dir           => $self->output_dir,
					-options               => $self->param('options'),
					-program               => $self->param('program_file'),
					-job_name              => $self->job_name,
					-samtools_path         => $samtools_path,
                                        -output_format         => $output_format,
                                        -bedGraphToBigWig_path => $bedGraph_to_bigWig_path,
                                        -chrom_sizes_file      => $chrom_sizes_file,
                                        -dedupe                => $dedup, 
                                        -mcr_root              => $mcr_root,
                              );


  $self->run_program( $wiggler ); 
  $self->output_param( 'bedgraph', ${$wiggler->output_files}[0] ) if $output_format eq 'bg'; 
  $self->output_param( 'bigwig', ${$wiggler->output_files}[0] )   if $output_format eq 'bw'; 
}


1;
