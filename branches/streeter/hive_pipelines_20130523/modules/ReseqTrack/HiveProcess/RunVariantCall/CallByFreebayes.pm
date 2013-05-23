
package ReseqTrack::HiveProcess::RunVariantCall::CallByFreebayes;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::RunVariantCall::CallByFreebayes;
use ReseqTrack::HiveUtils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $bam = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $fai = $self->param('fai') || throw("'fai' is an obligatory parameter");
    my $bed = $self->param('bed');
    my $SQ_start = $self->param('SQ_start');
    my $SQ_end = $self->param('SQ_end');
    my $bp_start = $self->param('bp_start');
    my $bp_end = $self->param('bp_end');
    my $overlap = $self->param('overlap') || 0;
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
      $slices = bed_to_fai(bed => $bed, parent_slices => $slices);
    }
    throw("expected exactly one slice, not ". scalar @$slices) if scalar @$slices != 1;
    $slices->[0]->extend($overlap);

    $self->data_dbc->disconnect_when_inactive(1);

    my $variant_caller = ReseqTrack::Tools::RunVariantCall::CallByFreebayes->new(
          -input_files => $bam,
          -program => $self->param('freebayes'),
          -bgzip => $self->param('bgzip'),
          -working_dir => $self->output_dir,
          -job_name => $self->job_name,
          -chrom => $slices->[0]->SQ_name,
          -reference => $self->param('reference'),
          -options => $self->param('options'),
          );
    if (!$slices->[0]->is_whole_SQ) {
      $variant_caller->region_start($slices->[0]->start);
      $variant_caller->region_end($slices->[0]->end);
    }
    $variant_caller->run;

    $self->output_this_branch('vcf' => $variant_caller->output_files);
    $self->data_dbc->disconnect_when_inactive(0);

}

1;

