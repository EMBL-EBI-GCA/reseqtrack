package ReseqTrack::Tools::RunVariantCall::CallByFreebayes;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);

use base qw(ReseqTrack::Tools::RunVariantCall);

sub DEFAULT_OPTIONS { return {
        };
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $bgzip)
    = rearrange([ qw(
    BGZIP
    )], @args);    
  
  ## Set defaults
  $self->program('freebayes') if (! $self->program);
  $self->bgzip($bgzip|| 'bgzip');

  return $self;
}

sub run_program {
    my ($self) = @_;

    throw("do not have have a reference") if !$self->reference;
    check_file_exists($self->reference);
    check_executable($self->bgzip);

    my $output_vcf = $self->working_dir .'/'. $self->job_name . '.vcf.gz';
    $output_vcf =~ s{//}{/};
    
    my @cmd_words;
    push(@cmd_words, $self->program);
    if ( $self->options ) {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;
        push(@cmd_words, "--$tag", $value);
      }
    }        
    push(@cmd_words, '-f', $self->reference);

    if (my $region_string = $self->chrom) {
      if (defined $self->region_start && defined $self->region_end) {
        $region_string .= ':' . $self->region_start . '..' . $self->region_end;
      }
      push(@cmd_words, '-r', $region_string);
    }
    foreach my $bam (@{$self->input_files}) {
      push(@cmd_words, '-b', $bam);
    }

    push(@cmd_words, '|', $self->bgzip, '-c');
    push(@cmd_words, '>', $output_vcf);

    my $cmd = join(' ', @cmd_words);
    $self->output_files($output_vcf);
    $self->execute_command_line ($cmd);

    return $self;

}


sub bgzip{
  my ($self, $bgzip) = @_;
  if ($bgzip) {
    $self->{'bgzip'} = $bgzip;
  }
  return $self->{'bgzip'};
}

1;

