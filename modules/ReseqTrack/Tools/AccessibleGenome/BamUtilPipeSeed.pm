package AccessibleGenome::BamUtilPipeSeed;

use ReseqTrack::Tools::Exception;
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

	my @tmp = split(/\./, $basename);
	my $pop = $tmp[3];
	my $sample = $tmp[0];
	
	$output_hash->{POP} = $pop;
	$output_hash{SAMPLE} = $sample;
	$output_hash->{fai} = '/nfs/production/reseq-info/work/reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai';
	
	push(@seed_params, $seed_params);

    $output_hash->{bam} = $name;

    # Temporary line:
    #last SEED if @seed_params >=5;
    #last SEED;
  }

  $self->seed_params(\@seed_params);
}

1;
