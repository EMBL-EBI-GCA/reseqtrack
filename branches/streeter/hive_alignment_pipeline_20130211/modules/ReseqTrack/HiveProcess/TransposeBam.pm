
package ReseqTrack::HiveProcess::TransposeBam;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::RunTransposeBam;
use ReseqTrack::HiveUtils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $bams = $self->param('bams') || throw("'bam' is an obligatory parameter");
    my $fai = $self->param('fai') || throw("'fai' is an obligatory parameter");
    my $bed = $self->param('bed');
    my $SQ_start = $self->param('SQ_start');
    my $SQ_end = $self->param('SQ_end');
    my $bp_start = $self->param('bp_start');
    my $bp_end = $self->param('bp_end');
    my $overlap = $self->param('region_overlap') || 0;
    throw("'SQ_start' is an obligatory parameter") if !defined $SQ_start;
    throw("'SQ_end' is an obligatory parameter") if !defined $SQ_end;
    throw("'bp_start' is an obligatory parameter") if !defined $bp_start;
    throw("'bp_end' is an obligatory parameter") if !defined $bp_end;

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

    $self->data_dbc->disconnect_when_inactive(1);

    my $bam_transposer = ReseqTrack::Tools::RunTransposeBam->new(
          -input_files => $bams,
          -program => $self->param('program_file'),
          -working_dir => $self->output_dir,
          -job_name => $self->job_name,
          -regions => \@regions,
          -options => {'build_index' => $self->param('create_index')},
          );
    $bam_transposer->run;

    my %output_data;
    $output_data{'bam'} = [grep { /\.bam$/ } @{$bam_transposer->output_files}];
    if ($self->param('create_index')) {
      $output_data{'bai'} = [grep { /\.bai$/ } @{$bam_transposer->output_files}];
    }

    $self->output_this_branch(%output_data);
    $self->data_dbc->disconnect_when_inactive(0);

}

1;

