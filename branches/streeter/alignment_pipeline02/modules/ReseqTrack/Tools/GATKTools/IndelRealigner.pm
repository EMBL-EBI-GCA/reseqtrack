=pod

=head1 NAME

ReseqTrack::Tools::GATKTools::IndelRealigner

=head1 SYNOPSIS

Object to create a bam file that in realigned around
known indel sites using GATK RealignerTargetCreator and
IndelRealigner

example

my $REALIGN_AROUND_INDELS = $IR->new(
		     -reference       => $reference,
		     -input_files     => $input{input_files},
		     -rtc_knowns      =>" -known $millsindels.vcf.gz ",
		     -working_dir     => $input{working_dir},
);



=cut

package ReseqTrack::Tools::GATKTools::IndelRealigner;

use strict;
use warnings;

use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);

use base qw(ReseqTrack::Tools::GATKTools);

sub DEFAULT_OPTIONS { return {
        'threads' => 1,
        'lod' => 0.4,
        'model' => 'KNOWNS_ONLY',
        };
}

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

        #setting defaults
        if (!$self->jar_file) {
          $self->jar_file ("GenomeAnalysisTK.jar");
        }

	return $self;
}


sub run_program {
	my $self = shift;

        throw "no input bam" if (!$self->input_bam);
        warn "No known indels files"
            if (! @{$self->known_sites_files});
        check_file_exists($_) foreach (@{$self->known_sites_files});
        $self->check_jar_file_exists;
        check_file_exists($self->reference);
        $self->check_bai_exists();


	$self->create_target_intervals_file();
	$self->create_indel_realign_bam();

	return;
}

sub create_target_intervals_file {
  my $self = shift;

  my $interval_file = $self->working_dir . '/'
        . $self->job_name . '.bam.interval_list';
  $interval_file =~ s{//}{/}g;

  my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
  push(@cmd_words, $self->gatk_path . '/' . $self->jar_file);
  push(@cmd_words, '-T', 'RealignerTargetCreator');
  push(@cmd_words, '-nt ' . $self->options('threads') || 1);

  foreach my $vcf (@{$self->known_sites_files}) {
    push(@cmd_words, '-known', $vcf);
  }
  push(@cmd_words, '-R', $self->reference);
  push(@cmd_words, '-I', $self->input_bam);
  push(@cmd_words, '-o', $interval_file);

  my $cmd = join(' ', @cmd_words);

  $self->intervals_file($interval_file);
  $self->execute_command_line ($cmd);

  return;
}


sub create_indel_realign_bam {
  my $self = shift;

  my $realigned_bam = $self->working_dir . '/'
        . $self->job_name. '.indel_realigned.bam';

  my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
  push(@cmd_words, $self->gatk_path . '/' . $self->jar_file);
  push(@cmd_words, '-T', 'IndelRealigner');
  push(@cmd_words, '-nt ' . $self->options('threads') || 1);
  push(@cmd_words, '-LOD', $self->options('lod')) if ($self->options('lod'));
  push(@cmd_words, '-model', $self->options('model')) if ($self->options('model'));
  push(@cmd_words, '--disable_bam_indexing');
  push(@cmd_words, '--targetIntervals', $self->intervals_file);

  foreach my $vcf (@{$self->known_sites_files}) {
    push(@cmd_words, '-known', $vcf);
  }
  push(@cmd_words, '-R', $self->reference);
  push(@cmd_words, '-I', $self->input_bam);
  push(@cmd_words, '-o', $realigned_bam);

  my $cmd = join(' ', @cmd_words);

  $self->output_files($realigned_bam);
  $self->execute_command_line ($cmd);

  return;

}

sub intervals_file {

  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{intervals_file} = $arg;
    $self->created_files($arg);
  }
  return $self->{intervals_file};
}

1;
