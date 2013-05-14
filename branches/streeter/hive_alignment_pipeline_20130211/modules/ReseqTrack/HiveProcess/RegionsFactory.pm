
package ReseqTrack::HiveProcess::RegionsFactory;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::HiveUtils::SequenceRegionsUtils qw(db_to_arrays fai_to_db);
use POSIX qw(ceil);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $command = $self->param('command') || throw("'command' is an obligatory parameter");

    if ($command eq 'load_db') {
      $self->load_db;
    }
    elsif ($command eq 'seq_region_factory') {
      $self->seq_region_factory;
    }
    elsif ($command eq 'window_factory') {
      $self->window_factory;
    }
    elsif ($command eq 'window_factory_with_targets') {
      $self->window_factory_with_targets;
    }
    else {
      throw("did not recognise command $command");
    }
}

#max_sequences overrules min_bases
sub seq_region_factory {
  my ($self) = @_;
  my $min_bases = $self->param('min_bases') || 0;
  my ($sequence_names, $sequence_lengths) = db_to_arrays($self->data_dbc);
  my $num_sequences = scalar @$sequence_names;
  throw("no sequences in database") if !$num_sequences;
  my $max_sequences = $self->param('max_sequences') || $num_sequences;
  my @regions;
  my $start_index = 0;
  my $base_counter = 0;
  my $end_index = 0;
  SEQ:
  while (1) {
    if ($end_index == $num_sequences-1) {
      push(@regions, [$start_index, $end_index]);
      last SEQ;
    }
    $base_counter += $sequence_lengths->[$end_index];
    if ($base_counter < $min_bases && $end_index - $start_index +1 < $max_sequences) {
      $end_index += 1;
      next SEQ;
    }
    push(@regions, [$start_index, $end_index]);
    $start_index = $end_index + 1;
    $end_index = $start_index;
    $base_counter = 0;
  }
  foreach my $region (@regions) {
    my $label = $sequence_names->[$region->[0]];
    if ($region->[0] != $region->[1]) {
      $label .= '-' . $sequence_names->[$region->[1]];
    }
    $self->output_child_branches('seq_index_start' => $region->[0], 'seq_index_end' => $region->[1], 'label' => $label);
  }
}

sub window_factory {
  my ($self) = @_;
  my $max_bases = $self->param('max_bases') || throw('max_bases is an obligatory parameter');
  my $seq_index_start = $self->param('seq_index_start') || 0;
  my $seq_index_end = $self->param('seq_index_end');
  my ($sequence_names, $sequence_lengths) = db_to_arrays($self->data_dbc);
  my $num_sequences = scalar @$sequence_names;
  throw("no sequences in database") if !$num_sequences;
  $seq_index_end = $num_sequences -1 if !defined $seq_index_end;
  foreach my $seq_index ($seq_index_start..$seq_index_end) {
    my $seq_length = $sequence_lengths->[$seq_index];
    my $seq_name = $sequence_names->[$seq_index];
    my $num_regions = ceil( $seq_length / $max_bases );
    my $region_length = ceil( $seq_length / $num_regions);
    foreach my $i (1..$num_regions) {
      my $region_end = $i * $region_length;
      my $region_start = $region_end - $region_length + 1;
      if ($i == $num_regions) {
        $region_end = $seq_length;
      }
      my $label = "$seq_name.$region_start-$region_end";
      $self->output_child_branches('seq_index' => $seq_index, 'region_start' => $region_start, 'region_end' => $region_end, 'label' => $label);
    }
  }

}

# similar to window_factory, but only creates window in regions allowed by the bed file
sub window_factory_with_targets {
  my ($self) = @_;
  my $max_bases = $self->param('max_bases') || throw('max_bases is an obligatory parameter');
  my $bed_file = $self->param('targets_bed_file') || throw('targets_bed_file is an obligatory parameter with the command window_factory_with_targets');
  my $seq_index_start = $self->param('seq_index_start') || 0;
  my $seq_index_end = $self->param('seq_index_end');
  my ($sequence_names, $sequence_lengths) = db_to_arrays($self->data_dbc);
  throw("no sequences in database") if ! scalar @$sequence_names;
  my %allowed_seq_name_index;
  foreach my $seq_index ($seq_index_start..$seq_index_end) {
    $allowed_seq_name_index{$sequence_names->[$seq_index]} = $seq_index;
  }

  open my $BED, '<', $bed_file or throw("could not open $bed_file $!");
  LINE:
  while (my $line = <$BED>) {
    chomp $line;
    my ($target_seq_name, $target_start, $target_end) = split("\t", $line);
    next LINE if !defined $allowed_seq_name_index{$target_seq_name};
    my $target_length = $target_end - $target_start +1;
    my $num_regions = ceil($target_length / $max_bases);
    my $region_length = ceil($target_length / $num_regions);
    foreach my $i (0..$num_regions-1) {
      my $region_start = $target_start + $i*$region_length;
      my $region_end = ($i == $num_regions -1) ? $target_end : $region_start + $region_length -1;
      my $label = "$target_seq_name.$region_start-$region_end";
      $self->output_child_branches('seq_index' => $allowed_seq_name_index{$target_seq_name}, 'region_start' => $region_start, 'region_end' => $region_end, 'label' => $label);
    }
  }
  close BED;

}

sub load_db {
  my ($self) = @_;
  my $fai = $self->param('fai') || throw("need a fai to load db");
  my $hive_dbc = $self->data_dbc();
  fai_to_db($hive_dbc, $fai);
}

1;

