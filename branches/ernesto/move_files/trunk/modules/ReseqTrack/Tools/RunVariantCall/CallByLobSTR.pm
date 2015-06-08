package ReseqTrack::Tools::RunVariantCall::CallByLobSTR;

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

  my ( $str_info, $noise_model, $ref_index_prefix, $bgzip)
    = rearrange([ qw(
    STR_INFO NOISE_MODEL REF_INDEX_PREFIX BGZIP
    )], @args);    
  
  $self->str_info($str_info);
  $self->noise_model($noise_model);
  $self->ref_index_prefix($ref_index_prefix);

  ## Set defaults
  $self->program('allelotype') if (! $self->program);
  $self->bgzip($bgzip|| 'bgzip');

  return $self;
}

sub run_program {
    my ($self) = @_;

    my $str_info = $self->str_info;
    my $noise_model = $self->noise_model;
    my $ref_index_prefix = $self->ref_index_prefix;

    throw("do not have have a ref_index_refix") if !$ref_index_prefix;
    throw("do not have have a noise model") if !$noise_model;
    check_executable($self->bgzip);
    check_file_exists($str_info);
    check_file_exists("$noise_model.stepmodel");
    check_file_exists("$noise_model.stuttermodel");
    check_file_exists("$noise_model.stutterproblem");

    if (defined $self->region_start || defined $self->region_end) {
      throw("LobSTR does not take a chromosome region as input");
    }

    my $output_prefix = $self->working_dir .'/'. $self->job_name;
    $output_prefix =~ s{//}{/};
    
    my @cmd_words;
    push(@cmd_words, $self->program, '--command', 'classify');
    if ( $self->options ) {
      OPTION:
      while (my ($tag, $value) = each %{$self->options}) {
        next OPTION if !defined $value;
        push(@cmd_words, "--$tag", $value);
      }
    }        
    push(@cmd_words, '--index-prefix', $ref_index_prefix);
    push(@cmd_words, '--noise_model', $noise_model);
    push(@cmd_words, '--strinfo', $str_info);

    if (my $chrom = $self->chrom) {
      push(@cmd_words, '--chrom', $chrom);
    }
    push(@cmd_words, '--bam', join(',', @{$self->input_files}));

    push(@cmd_words, '--out', $output_prefix);

    my $cmd = join(' ', @cmd_words);
    $self->created_files("$output_prefix.allelotype.stats");
    $self->output_files("$output_prefix.vcf");
    $self->execute_command_line($cmd);

    return $self;

}


sub ref_index_prefix{
  my ($self, $ref_index_prefix) = @_;
  if ($ref_index_prefix) {
    $self->{'ref_index_prefix'} = $ref_index_prefix;
  }
  return $self->{'ref_index_prefix'};
}

sub str_info{
  my ($self, $str_info) = @_;
  if ($str_info) {
    $self->{'str_info'} = $str_info;
  }
  return $self->{'str_info'};
}

sub bgzip{
  my ($self, $bgzip) = @_;
  if ($bgzip) {
    $self->{'bgzip'} = $bgzip;
  }
  return $self->{'bgzip'};
}


sub noise_model{
  my ($self, $noise_model) = @_;
  if ($noise_model) {
    $self->{'noise_model'} = $noise_model;
  }
  return $self->{'noise_model'};
}

1;

