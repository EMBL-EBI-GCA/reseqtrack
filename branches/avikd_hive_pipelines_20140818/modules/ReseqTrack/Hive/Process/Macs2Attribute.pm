package ReseqTrack::Hive::Process::Macs2Attribute;

use strict;
use warnings;
use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::GeneralUtils qw( execute_system_command execute_pipe_system_command );
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use Statistics::Descriptive;

sub run {
  my $self = shift @_;
  
  $self->param_required('bam'); 
  my $bed           = $self->param_required('bed'); 
  my $bams         = $self->param_as_array( 'bam' );
  my $db_param = $self->param_required('reseqtrack_db');
  my $samtools   = $self->param_required('samtools');
  my $bedtools    = $self->param_required('bedtools');
  
  throw('expecting single bed file') if ( ref($bed) eq 'ARRAY' );
  check_file_exists($bed);
  
  foreach my $bam (@$bams) {
    check_file_exists($bam);
  }
  my $bam = $$bams[0];
  
  my %peak_stats = _gather_peak_stats( $bed, $bam, $samtools, $bedtools );

  $self->output_param('attribute_metrics',  \%peak_stats );
}

sub _gather_peak_stats {
    my ( $bed_file, $bam_file, $samtools_path, $bedtools_path ) = @_;

    my %peak_stats;

    $peak_stats{total_reads} = _get_total_reads( $bam_file, $samtools_path );
    $peak_stats{reads_in_peaks} =
      _get_reads_in_peaks( $bam_file, $bed_file, $bedtools_path );

    if ( $peak_stats{total_reads} ) {
        $peak_stats{peak_enrichment} =
          ( $peak_stats{reads_in_peaks} / $peak_stats{total_reads} );
    }
    else {
        $peak_stats{peak_enrichment} = undef;
    }

    my $stats = Statistics::Descriptive::Full->new();

    my $bed_fh;
    if ( $bed_file =~ m/\.gz$/ ) {
        open $bed_fh, '-|', 'gzip', '-dc', $bed_file;
    }
    else {
        open( $bed_fh, '<', $bed_file );
    }
    while (<$bed_fh>) {
        next if (m/^#/);
        chomp;
        my ( $seq, $start, $end ) = split /\t/;
        my $length = $end - $start + 1;
        $stats->add_data($length);
    }
    close $bed_fh;

    $peak_stats{region_count}           = $stats->count();
    $peak_stats{region_length_mean}     = $stats->mean();
    $peak_stats{region_length_median}   = $stats->median();
    $peak_stats{region_length_std_dev}  = $stats->standard_deviation();
    $peak_stats{region_length_variance} = $stats->variance();
    $peak_stats{region_length_total}    = $stats->sum();

    return %peak_stats;
}

sub _get_total_reads {
    my ( $bam_file, $samtools_path ) = @_;

    my $read_count =
      execute_pipe_system_command("$samtools_path view -c $bam_file");

    return $read_count;
}

sub _get_reads_in_peaks {
    my ( $bam_file, $bed_file, $bedtools_path ) = @_;

    my $reads_in_peaks = execute_pipe_system_command(
        "$bedtools_path intersect -a $bam_file -b $bed_file -bed | wc -l");

    return $reads_in_peaks;
}

1;

