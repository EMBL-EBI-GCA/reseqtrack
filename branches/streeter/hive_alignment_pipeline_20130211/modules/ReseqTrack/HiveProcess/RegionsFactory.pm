
package ReseqTrack::HiveProcess::RegionsFactory;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use POSIX qw(floor);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $fai = $self->param('fai') || die "'fai' is an obligatory parameter";
    my $child_num_bases = $self->param('child_num_bases');
    my $child_num_regions = $self->param('child_num_regions') || 1;
    my $num_bases_tolerance = $self->get_param_array('num_bases_tolerance');
    my $whole_SQ = $self->param('whole_SQ');
    my $parent_region_strings = $self->get_param_array('parent_regions');

    if (!defined $num_bases_tolerance) {
      $num_bases_tolerance = floor(0.1 * $child_num_bases);
    }

    my @fai_data;
    my %SQ_fai_order;
    open my $FAI, '<', $fai or throw("could not open $fai: $!");
    LINE:
    while (my $line = <$FAI>) {
      my ($SQ, $length) = split(/\s+/, $line);
      throw("did not recognise line in fai: $line") if (!defined $length);
      push(@fai_data, [$SQ, $length]);
      $SQ_fai_order{$SQ} = $#fai_data;
    }
    close $FAI;

    #my @parent_regions = map {[ /([^:]+)(?::(\d+)-(\d+))?/ ]} @$parent_region_strings;
    my @sorted_parent_regions = sort {$SQ_fai_order{$a->[0]} <=> $SQ_fai_order{$b->[0]}
                                  || $a->[1] <=> $b->[1]} map {[ /([^:]+)(?::(\d+)-(\d+))?/ ]} @$parent_region_strings;


    my %parent_regions;
    foreach my $region_string (@$parent_region_strings) {
      my ($SQ, $start, $end) = $region_string =~ /([^:]+)(?::(\d+)-(\d+))?/;
      if (!defined $start) {
        $start = 1;
        $end = $fai_data[$SQ_fai_order{$SQ}]->[1];
      }
      push(@{$parent_regions{$SQ}}, [$start, $end]);
    }

    my $num_groups = 0;
    my $base_counter = 0;
    my @unprocessed_regions;
    SQ:
    foreach my $fai_SQ_data (@fai_data) {
      my ($SQ, $SQ_length) = @$fai_SQ_data;
      my ($current_start, $current_end) = (1,0);

      ALLOWED_REGION:
      while (1) {
        my $whole_SQ_is_allowed = 0;
        my $allowed_end;
        if (scalar @$parent_region_strings) {
          last SQ if !@sorted_parent_regions;
          next SQ if $sorted_parent_regions[0]->[0] ne $SQ;
          $current_start = $sorted_parent_regions[0]->[1];
          $allowed_end = $sorted_parent_regions[0]->[2];
          throw("have region end but not region start $SQ $current_start") if (defined $current_start && !defined $allowed_end);
        }
        if (!defined $allowed_end) {
          $whole_SQ_is_allowed = 1;
          $allowed_end = $SQ_length;
        }

        OUTPUT_REGION:
        while (1) {
          if ($current_end >= $SQ_length) {
            if (!$child_num_bases || $base_counter >= $child_num_bases || scalar @unprocessed_regions > $child_num_regions) {
              process_regions(\@unprocessed_regions, \$num_groups);
            }
            last OUTPUT_REGION;
          }
          if ($whole_SQ_is_allowed) {
            if ($whole_SQ || !$child_num_bases || $child_num_bases + $num_bases_tolerance >= $base_counter + $SQ_length) {
              push(@unprocessed_regions, [$SQ]);
              $base_counter += $SQ_length;
              $current_end = $SQ_length;
              redo OUTPUT_REGION;
            }
            else {
              $current_start = 1;
            }
          }

          if ($child_num_bases) {
            $current_end = $current_start + $child_num_bases - $base_counter - 1;
            $current_end = $current_end >= $allowed_end - $num_bases_tolerance ? $allowed_end : $current_end;
          }
          else {
            $current_end = $allowed_end;
          }
          push(@unprocessed_regions, [$SQ, $current_start, $current_end]);
        }

        last ALLOWED_REGION if !scalar @$parent_region_strings;
        shift @sorted_parent_regions;
      }

    }

}

sub process_regions {
  my ($self, $regions_group, $group_num_ref) = @_;

  my $label;
  if (scalar @$regions_group > 1) {
    $$group_num_ref += 1;
    $label = "SQgroup$$group_num_ref";
  }
  else {
    $label = join('_', map{defined $_} @{$regions_group->[0]}{qw(SQ start end)});
  }
  my @regions_strings;
  foreach my $region ($regions_group) {
    my ($SQ, $start, $end) = @$region;
    my $string = $SQ;
    $string .= ":$start" if defined $start;
    $string .= "-$end" if defined $end;
    push(@regions_strings, $string);
  }

  $self->output_child_branches('region' => \@regions_strings, 'label' => $label);
}

1;

