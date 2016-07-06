package AccessibleGenome::PipeSeed;
use base ('ReseqTrack::Hive::PipeSeed::BasePipeSeed');

use File::Basename qw(fileparse);

sub create_seed_params {
  my ($self) = @_;
  ReseqTrack::Hive::PipeSeed::BasePipeSeed::create_seed_params(@_);

  my @seed_params;

  SEED:
  foreach my $seed_params (@{$self->seed_params}) {
    my ($file, $output_hash) = @$seed_params;

    my $name = $file->name;
    my ($basename) = fileparse($name, '.bam');
    $output_hash->{file_basename} = $basename;

    # Temporary line:
    #next SEED if $basename !~ /^IRBT/;

    my ($nextgen_group) = $name =~ /\/([A-Z]{4})-[A-Z]+\d+-\d+\//;
    #$nextgen_group = 'AUCH';
    $output_hash->{NEXTGEN_GROUP} = $nextgen_group;

    if ($basename =~ /OARv3_1/ ) {
      $output_hash->{fai} = '/nfs/1000g-work/G1K/work/streeter/nextgen/alignments/sheep_reference/OARv3.1/sheep_oarv3.1.20130206/sheep_oarv3.1.20130206.fa.fai';
      push(@seed_params, $seed_params);
    }
    elsif ($basename =~ /CHIR1_0/ ) {
      $output_hash->{fai} = '/nfs/1000g-work/G1K/work/streeter/nextgen/alignments/goat_reference/CHIR1.0/seq_accessions_20130222/goat_chir1.0.20130222.fa.fai';
      push(@seed_params, $seed_params);
    }
    elsif ($basename =~ /UMD3_1/ ) {
      $output_hash->{fai} = '/nfs/1000g-work/G1K/work/streeter/nextgen/alignments/cow_reference/UMD3.1/cow_umd3.1.20131028.fa.fai';
      push(@seed_params, $seed_params);
    }
    else {
      print "did not recognise basename $basename\n";
      next SEED;
    }

    $output_hash->{bam} = $name;

    # Temporary line:
    #last SEED if @seed_params >=5;
    #last SEED;
  }

  $self->seed_params(\@seed_params);
}

1;
