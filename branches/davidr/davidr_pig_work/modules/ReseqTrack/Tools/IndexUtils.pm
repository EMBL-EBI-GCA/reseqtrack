package ReseqTrack::Tools::IndexUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);


use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(get_md5hash_from_sequence_index
             get_md5hash_from_alignment_index
            );





sub get_md5hash_from_sequence_index{
  my ($index, $skip_withdrawn) = @_;
  my $index_hash = get_hash_from_tab_file($index, 0);
  my %md5_hash;
  foreach my $key(keys(%$index_hash)){
    my $string = $index_hash->{$key};
    my @values = split /\t/, $string;
    next if($values[20] && $skip_withdrawn);
    $md5_hash{$values[0]} = $values[1];
  }
  return \%md5_hash;
}

sub get_md5hash_from_alignment_index{
  my ($index) = @_;
  my $index_hash = get_hash_from_tab_file($index, 0);
  my %md5_hash;
  foreach my $key(keys(%$index_hash)){
    my $string = $index_hash->{$key};
    my @values = split /\t/, $string;
    $md5_hash{$values[0]} = $values[1];
    $md5_hash{$values[4]} = $values[5];
  }
  return \%md5_hash;
}






1;
