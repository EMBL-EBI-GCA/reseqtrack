
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
    else {
      throw("did not recognise command $command");
    }
}

sub seq_region_factory {
  my ($self) = @_;
  $min_bases = $self->param('min_bases') || 0;
  my ($sequence_names, $sequence_lengths) = $self->_db_to_arrays;
  my $num_sequences = scalar @$sequence_names;
  throw("no sequences in database") if !$num_sequences;
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
    $base_counter += $sequence_lengths[$end_index];
    if ($base_counter < $min_bases) {
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
  $max_bases = $self->param('max_bases') || throw('max_bases is an obligatory parameter');
  $seq_index_start = $self->param('seq_index_start') || 0;
  $seq_index_end = $self->param('seq_index_end');
  my ($sequence_names, $sequence_lengths) = $self->_db_to_arrays;
  my $num_sequences = scalar @$sequence_names;
  throw("no sequences in database") if !$num_sequences;
  $seq_index_end = $num_sequences -1 if !defined $seq_index_end;
  foreach my $seq_index ($seq_index_start..$seq_index_end) {
    my $seq_length = $sequence_lengths->[$seq_index];
    my $seq_name = $sequence_names->[$seq_index];
    my $num_regions = ceil( $seq_length / $max_bases );
    my $region_length = ceil( $seq_length / $num_regions);
    foreach my $i (1..$num_regions) {
      my %region = {'seq_index' => $seq_index, 'seq_name' => $seq_name};
      my $seq_end = $i * $region_length;
      my $seq_start = $region{'end'} - $region_length + 1;
      my $label = "$seq_name.$seq_start-$seq_end";
      $self->output_child_branches('seq_index' => $seq_index, 'seq_start' => $seq_start, 'seq_end' => $seq_end, 'label' => $label);
    }
  }

}

sub load_db {
  my ($self) = @_;
  my $fai = $self->param('fai') || throw("need a fai to load db");

  my $hive_dbc = $self->data_dbc();
  my $sql = "INSERT INTO sequence_region (sequence_region_index, name, length VALUES (?, ?, ?)";
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;

  my $sequence_region_index = 0;
  open my $FAI, '<', $fai or throw("could not open $fai: $!");
  while (my $line = <$FAI>) {
    my ($SQ, $length) = split(/\s+/, $line);
    throw("did not recognise line in fai: $line") if (!defined $length);
    $sth->bind_param(1, $sequence_region_index);
    $sth->bind_param(2, $SQ);
    $sth->bind_param(3, $length);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
  }
  close $FAI;
}

sub _db_to_arrays {
  my ($self) = @_;
  my $hive_dbc = $self->data_dbc();
  my @sequence_names;
  my @sequence_lengths;
  my $sql = "SELECT sequence_region_index, name, length FROM sequence_region";
  my $sth = $hive_dbc->prepare($sql) or die "could not prepare $sql: ".$hive_dbc->errstr;
  foreach my $row (@{$sth->fetchall_arrayref}) {
    $sequence_names[$row->[0]] = $row->[1];
    $sequence_lengths[$row->[0]] = $row->[2];
  }
  return \@sequence_names, \@sequence_lengths;
}

1;

