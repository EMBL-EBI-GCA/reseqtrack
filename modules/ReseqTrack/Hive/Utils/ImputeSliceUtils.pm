package ReseqTrack::Hive::Utils::ImputeSliceUtils;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use base qw(Exporter);

our @ISA = qw(Exporter);
our @EXPORT_OK = qw( reference_config_to_bed_slices get_chrX filter_haps_missing);

sub reference_config_to_bed_slices {
  my (%args) = @_;
  
  my ($reference_config, $SQ_start, $bp_start, $bp_end, $max_base, $gunzip, $chrX_string, $chrX_PAR1_start, $chrX_PAR1_end, $chrX_PAR2_start, $chrX_PAR2_end) = 
    @args{qw(reference_config SQ_start bp_start bp_end max_base gunzip chrX_string chrX_PAR1_start chrX_PAR1_end chrX_PAR2_start chrX_PAR2_end)};
    
  check_file_exists($reference_config);
  throw("CHROM name required") if(!$SQ_start);
  
  my $chrom_x = 0;
  
  if(defined $chrX_string && $SQ_start eq $chrX_string) {
    $chrom_x = 1;
    
    $SQ_start = get_chrX(
      SQ_start => $SQ_start, 
      bp_start => $bp_start, 
      bp_end => $bp_end, 
      chrX_PAR1_start => $chrX_PAR1_start, 
      chrX_PAR1_end => $chrX_PAR1_end, 
      chrX_PAR2_start => $chrX_PAR2_start, 
      chrX_PAR2_end => $chrX_PAR2_end,
    );
  }

  my @bed_slices;
  
  my $found_start = 0;
  
  open my $REF_CONF, '<', $reference_config or throw("could not open $reference_config: $!");
  LINE:
  while (my $line = <$REF_CONF>) {
    next LINE if($line=~ /^#/);
    
    my ($SQ, $map, $haps, $legends) = split(/\s+/, $line);
    
    throw("did not recognise line in reference_config: $line") if (!defined $legends);
    
    if($SQ eq $SQ_start){
      $found_start = 1;
      
      check_file_exists($legends);
      
      $SQ_start = $chrX_string if($chrom_x == 1);
      
      my $bed_slice = get_legends_slice($SQ, $legends, $max_base, $gunzip);
      push(@bed_slices, $bed_slice);
    }
    last LINE if ($found_start == 1);
  }
  close $REF_CONF;
  
  throw("CHROM $SQ_start not found in reference_config file") if($found_start == 0);
  
  return \@bed_slices;
}

sub get_legends_slice {
  my ($SQ, $legends, $max_base, $gunzip) = @_;
  
  my @slice;
  my $LEGENDS;
  
  if($legends =~ /\.gz$/) {
    open $LEGENDS, "$gunzip -cd $legends |" or throw("could not open $legends: $!");
  }
  else  {
    open $LEGENDS, '<', $legends or throw("could not open $legends: $!");
  }
  
  my $count = 0;
  my $start_pos;
  my $end_pos;
      
  my $line=<$LEGENDS>; ## removing header
  
  while( $line = <$LEGENDS> ) {
    chomp($line);
    
    my ($rs_id, $pos) = split(/\s+/, $line);

    throw("did not recognise line in reference_config: $line") if (!defined $pos);
    
    $count++;
    
    if($count == 1)
    {
      $start_pos = $pos;
      $end_pos = $pos;
    }
    my $distance = $end_pos - $start_pos;
    
    if( $distance > $max_base ) {
      $count = 1;
      push @slice, $SQ ." ". $start_pos . " " . $end_pos;
      $start_pos = $pos;
      $end_pos = $pos;
    }
    else {
      $end_pos = $pos;
    } 
  }
  close $LEGENDS;
  push @slice, $SQ ." ". $start_pos . " " . $end_pos;
  
  return \@slice;
}

sub get_chrX {
    my (%args) = @_;
    
    my ($SQ_start, $bp_start, $bp_end, $chrX_PAR1_start, $chrX_PAR1_end, $chrX_PAR2_start, $chrX_PAR2_end ) = 
        @args{qw( SQ_start bp_start bp_end chrX_PAR1_start chrX_PAR1_end chrX_PAR2_start chrX_PAR2_end )};
  
  if($bp_start >= $chrX_PAR1_start && $bp_end <= $chrX_PAR1_end ){
          $SQ_start = "X_PAR1";
        }
  elsif($bp_start > $chrX_PAR1_end && $bp_end < $chrX_PAR2_start){
          $SQ_start = "X_nonPAR";
        }
  elsif($bp_start >= $chrX_PAR2_start && $bp_end <= $chrX_PAR2_end){
          $SQ_start = "X_PAR2";
        }
  return $SQ_start;      
    
}

sub filter_haps_missing {
  my (%args) = @_;
  
  my($haps, $slice, $SQ_start) =  @args{qw(haps slice SQ_start)};
  
  my @new_slice;
  
  my %slice_hash;
  
  foreach my $slice_line(@$slice) {
  
    my($chrom, $start, $end) = split(/\s+/, $slice_line);
    
    throw("Chrom mismatch: slice $chrom and haps  $SQ_start") unless($SQ_start eq $chrom);
    
    $slice_hash{$start." ".$end} = 0;
  }
  
  my $HAPS;
  
  open $HAPS, '<', $haps or throw("could not open $haps: $!");
  
  while( my $line = <$HAPS> ) {
    chomp($line);
    
    my ($id_1, $id_2, $pos) = split(/\s+/, $line);

    throw("did not recognise line in $haps : $line") if (!defined $pos);
    
    foreach my $slice_key (keys %slice_hash) {
      my ($start_pos, $end_pos) = split(/\s+/, $slice_key);
      
      $slice_hash{$slice_key}++ if($pos >= $start_pos && $pos <= $end_pos);   
    }
  }
  close $HAPS;
  
  foreach my $slice_key (keys %slice_hash) {
    if($slice_hash{$slice_key} > 0 ) {
      
      push @new_slice, $SQ_start." ".$slice_key;
    }
  }
  return \@new_slice;
}
1;