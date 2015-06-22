
package ReseqTrack::Hive::Utils::SequenceSliceUtils;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Hive::Utils::SequenceSlice;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use base qw(Exporter);

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(fai_to_slices bed_to_slices);

sub fai_to_slices {
  my (%args) = @_;
  my ($fai, $SQ_start, $SQ_end, $bp_start, $bp_end) = 
    @args{qw(fai SQ_start SQ_end bp_start bp_end)};
  check_file_exists($fai);

  my @slices;
  my $found_start = defined $SQ_start ? 0 : 1;
  open my $FAI, '<', $fai or throw("could not open $fai: $!");
  LINE:
  while (my $line = <$FAI>) {
    my ($SQ, $length) = split(/\s+/, $line);
    throw("did not recognise line in fai: $line") if (!defined $length);
    if (!$found_start && $SQ eq $SQ_start) {
      $found_start = 1;
    }
    next LINE if !$found_start;

    my $slice = ReseqTrack::Hive::Utils::SequenceSlice->new(
        -SQ_name => $SQ,
        -SQ_length => $length,
        );
    push(@slices, $slice);

    last LINE if (defined $SQ_end && $SQ eq $SQ_end);
  }
  close $FAI;

  if (scalar @slices && defined $bp_start) {
    $slices[0]->start($bp_start);
  }
  if (scalar @slices && defined $bp_end) {
    $slices[-1]->end($bp_end);
  }

  return \@slices;
}

sub bed_to_slices {
  my (%args) = @_;
  my ($bed, $parent_slices) = 
    @args{qw(bed parent_slices)};
  throw("no parent slices") if ! @$parent_slices;
  check_file_exists($bed);
  
  my $num_slices = scalar @$parent_slices;

  my @bed_slices;
  open my $BED, '<', $bed or throw("could not open $bed: $!");
  LINE:
  while (my $line = <$BED>) {
    chomp($line);
    
    
    my ($SQ, $start, $end) = split(/\s+/, $line);
    throw("did not recognise line in bed: $line") if (!defined $end);
    my @parent_index = grep {$parent_slices->[$_]->end >= $start}
                       grep {$parent_slices->[$_]->start <= $end}
                       grep {$parent_slices->[$_]->SQ_name eq $SQ} (0..$num_slices-1);
    next LINE if !@parent_index;
    throw("found multiple parent slices") if scalar @parent_index !=1;
    my $parent_slice = $parent_slices->[$parent_index[0]];
    if ($start < $parent_slice->start) {
      $start = $parent_slice->start;
    }
    if ($end > $parent_slice->end) {
      $end = $parent_slice->end;
    }
    
    my $bed_slice = ReseqTrack::Hive::Utils::SequenceSlice->new(
        -SQ_name => $SQ,
        -SQ_length => $parent_slice->SQ_length,
        -start => $start,
        -end => $end,
        );
    push(@{$bed_slices[$parent_index[0]]}, $bed_slice);
   
  }
  close $BED;
  
  
  my @sorted_bed_slices = map {sort {$a->start <=> $b->start} @$_} @bed_slices;
    return \@sorted_bed_slices;
}


1;
