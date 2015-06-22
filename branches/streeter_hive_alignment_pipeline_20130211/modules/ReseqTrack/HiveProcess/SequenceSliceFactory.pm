
package ReseqTrack::HiveProcess::SequenceSliceFactory;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::HiveUtils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);
use POSIX qw(ceil);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
  my ($self) = @_;
  my $fai = $self->param('fai') || throw("need a fai to to run SequenceSliceFactory");
  my $bed = $self->param('bed');
  my $max_sequences = $self->param('max_sequences') || 0;
  my $num_bases = $self->param('num_bases') || 0;
  my $SQ_start = $self->param('SQ_start');
  my $SQ_end = $self->param('SQ_end');
  my $bp_start = $self->param('bp_start');
  my $bp_end = $self->param('bp_end');

  my $slices = fai_to_slices(
          fai => $fai,
          SQ_start => $SQ_start, SQ_end => $SQ_end,
          bp_start => $bp_start, bp_end => $bp_end,
          );

  if (defined $bed) {
    $slices = bed_to_slices(bed => $bed, parent_slices => $slices);
  }

  if ($num_bases) {
    my @split_slices = map {@{$_->split($num_bases)}} @$slices;
    $slices = \@split_slices;
  }

  my $base_counter = 0;
  my $sequence_counter = 0;
  my @child_slices;
  SLICE:
  foreach my $slice (@$slices) {
    if ($base_counter == 0 || $base_counter + $slice->length > $num_bases
        || ($max_sequences && $sequence_counter >= $max_sequences) ) {
      push(@child_slices, [$slice, $slice]);
      $base_counter += $slice->length;
      $sequence_counter = 1;
    }
    else {
      $child_slices[-1][1] = $slice;
      $base_counter = $slice->length;
      $sequence_counter += 1;
    }
  }

  foreach my $child (@child_slices) {
    my $label;
    my $SQ_start = $child->[0]->SQ_name;
    my $SQ_end = $child->[1]->SQ_name;
    my $bp_start = $child->[0]->start;
    my $bp_end = $child->[1]->end;
    if ($SQ_start eq $SQ_end) {
      $label = $SQ_start;
      if ($bp_start != 1 || $bp_end != $child->[1]->SQ_length) {
        $label .= ".$bp_start-$bp_end";
      }
    }
    else {
      $label = $SQ_start;
      $label .= ".$bp_start" if $bp_start !=1;
      $label .= "-$SQ_end";
      $label .= ".$bp_end" if $bp_end != $child->[1]->SQ_length;
    }

    $self->output_child_branches(
          'SQ_start' => $child->[0]->SQ_name,
          'bp_start' => $child->[0]->start,
          'SQ_end' => $child->[1]->SQ_name,
          'bp_end' => $child->[1]->end,
          'label' => $label);
  }
}

1;

