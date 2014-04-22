
package ReseqTrack::Hive::Process::TransposeBam;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::RunTransposeBam;
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);
use ReseqTrack::Tools::Exception qw(throw);

sub param_defaults {
  return {
    bed => undef,
    region_overlap => 0,
    create_index => 1,
    uniquify_rg => 0,
  };
}


sub run {
    my $self = shift @_;

    $self->param_required('bam');
    my $bams = $self->param_as_array('bam');
    my $fai = $self->param_required('fai');
    my $bed = $self->param('bed');
    my $SQ_start = $self->param_required('SQ_start');
    my $SQ_end = $self->param_required('SQ_end');
    my $bp_start = $self->param_required('bp_start');
    my $bp_end = $self->param_required('bp_end');
    my $overlap = $self->param('region_overlap');
    my $create_index = $self->param('create_index');

    my $slices = fai_to_slices(
          fai => $fai,
          SQ_start => $SQ_start, SQ_end => $SQ_end,
          bp_start => $bp_start, bp_end => $bp_end,
          );
    if (defined $bed) {
      $slices = bed_to_slices(bed => $bed, parent_slices => $slices);
    }

    foreach my $slice (@$slices) {
      $slice->extend($overlap);
    }
    #$slices = join_overlapping_slices(slices => $slices, separation => 500);

    my @regions;
    foreach my $slice (@$slices) {
      my $region = $slice->SQ_name;
      if (!$slice->is_whole_SQ) {
        $region .= ':' . $slice->start . '-' . $slice->end;
      }
      push(@regions, $region);
    }

    my $bam_transposer = ReseqTrack::Tools::RunTransposeBam->new(
          -input_files => $bams,
          -program => $self->param('program_file'),
          -working_dir => $self->output_dir,
          -job_name => $self->job_name,
          -regions => \@regions,
          -options => {uniquify_rg => $self->param('uniquify_rg'), build_index => $self->param('create_index')},
          );

    $self->run_program($bam_transposer);

    $self->output_param('bam', $bam_transposer->output_bam_file);
    if ($create_index) {
      $self->output_param('bai', $bam_transposer->output_bai_file);
    }

}

1;

