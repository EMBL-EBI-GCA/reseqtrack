
package ReseqTrack::HiveUtils::SequenceRegionsUtils;

use strict;
use ReseqTrack::Tools::Exception qw(throw);
use base qw(Exporter);

our @ISA = qw(Exporter);
our @EXPORT_OK = qw(fai_to_db db_to_arrays get_region_strings get_sequence_name);

sub fai_to_db {
  my ($dbc, $fai) = @_;

  my $sql = "INSERT INTO sequence_region (sequence_region_index, name, length) VALUES (?, ?, ?)";
  my $sth = $dbc->prepare($sql) or die "could not prepare $sql: ".$dbc->errstr;

  my $sequence_region_index = 0;
  open my $FAI, '<', $fai or throw("could not open $fai: $!");
  while (my $line = <$FAI>) {
    my ($SQ, $length) = split(/\s+/, $line);
    throw("did not recognise line in fai: $line") if (!defined $length);
    $sth->bind_param(1, $sequence_region_index);
    $sth->bind_param(2, $SQ);
    $sth->bind_param(3, $length);
    $sth->execute() or die "could not execute $sql: ".$sth->errstr;
    $sequence_region_index +=1;
  }
  close $FAI;
}

sub db_to_arrays {
  my ($dbc) = @_;
  my @sequence_names;
  my @sequence_lengths;
  my $sql = "SELECT sequence_region_index, name, length FROM sequence_region";
  my $sth = $dbc->prepare($sql) or die "could not prepare $sql: ".$dbc->errstr;
  $sth->execute() or die "could not execute $sql: ".$sth->errstr;
  foreach my $row (@{$sth->fetchall_arrayref}) {
    $sequence_names[$row->[0]] = $row->[1];
    $sequence_lengths[$row->[0]] = $row->[2];
  }
  return \@sequence_names, \@sequence_lengths;
}

sub get_region_strings {
  my ($dbc, $seq_index_start, $seq_index_end, $position_start, $position_end) = @_;

  my ($sequence_names, $sequence_lengths) = db_to_arrays($dbc);

  my @region_strings;
  if ($seq_index_start == $seq_index_end) {
    my $string = $sequence_names->[$seq_index_start];
    if (defined $position_start || defined $position_end) {
      $string .= ':' . $position_start || 1;
      $string .=  '-' . defined $position_end ? $position_end : $sequence_lengths->[$seq_index_end];
    }
    push(@region_strings, $string);
  }
  else {
    my $first_string = $sequence_names->[$seq_index_start];
    if (defined $position_start) {
      $first_string .= ":$position_start-" . $sequence_lengths->[$seq_index_start];
    }
    push(@region_strings, $first_string);
    foreach my $seq_index ($seq_index_start+1..$seq_index_end-1) {
      push(@region_strings, $sequence_names->[$seq_index]);
    }
    my $last_string = $sequence_names->[$seq_index_end];
    if (defined $position_end) {
      $last_string .= ":1-$position_end";
    }
    push(@region_strings, $last_string);
  }
  return \@region_strings
}

sub get_sequence_name {
  my ($dbc, $seq_index) = @_;

  my ($sequence_names) = db_to_arrays($dbc);
  return $sequence_names->[$seq_index];
}

1;
