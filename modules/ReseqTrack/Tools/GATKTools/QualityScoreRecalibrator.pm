=pod

=head1 NAME

ReseqTrack::Tools::GATKTools::QualityScoreRecalibrator

=head1 SYNOPSIS

Object to create a bam file that in realigned around
known indel sites using GATK RealignerTargetCreator and
IndelRealigner

example

my $RECALIBRATOR = $QSR->new(
     -input_files     => '/path/to/bam'
     -java_exe        =>"/usr/bin/java" ,
     -jvm_args        =>"-Xmx4g",
     -GATK_PATH       =>$GATK/GenomeAnalysisTK/",
     -working_dir     => '/path/to/dir/',
     -reference       => '/path/to/reference',
     -known_sites_files => '/path/to/file',
);

=cut

package ReseqTrack::Tools::GATKTools::QualityScoreRecalibrator;

use strict;
use warnings;

use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);

use base qw(ReseqTrack::Tools::GATKTools);

sub DEFAULT_OPTIONS { return {
        'threads' => 1,
        'log_level' => 'INFO',
        'read_group_covariate' => 1,
        'quality_score_covariate' => 1,
        'cycle_covariate' => 1,
        'context_covariate' => 1,
	'intervals' => join(',', 1..22, 'X', 'Y', 'MT'),
        'solid_nocall_strategy' => 0,
        };
}

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	return $self;
}


sub run_program {
	my $self = shift;

        throw "no input bam" if (!$self->input_bam);
        check_file_exists($_) foreach (@{$self->known_sites_files});
        $self->check_jar_file_exists;
        check_file_exists($self->reference);
        $self->check_bai_exists();
        check_executable($self->java_exe);


	$self->create_recalibration_table();
	$self->create_recalibrated_bam();

	return;
}


sub create_recalibration_table {
  my $self = shift;

  my $covariates_file = $self->working_dir . '/'
        . $self->job_name . '.bam.recal_data.csv';
  $covariates_file =~ s{//}{/}g;

  my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
  push(@cmd_words, $self->gatk_path . '/' . $self->jar_file);
  push(@cmd_words, '-T', 'BaseRecalibrator');
  push(@cmd_words, '-nt ' . $self->options('threads') || 1);
  push(@cmd_words, '-l', $self->options('log_level')) if ($self->options('log_level'));
  push(@cmd_words, '-cov ReadGroupCovariate') if ($self->options('read_group_covariate'));
  push(@cmd_words, '-cov QualityScoreCovariate') if ($self->options('quality_score_covariate'));
  push(@cmd_words, '-cov CycleCovariate') if ($self->options('cycle_covariate'));
  push(@cmd_words, '-cov ContextCovariate') if ($self->options('context_covariate'));
  push(@cmd_words, '-allowPotentiallyMisencodedQuals') if ($self->options('allowPotentiallyMisencodedQuals'));

  if ($self->options('intervals')) {
    foreach my $interval (split(/,/, $self->options('intervals'))) {
	push(@cmd_words, '-L', $interval);
    }
  }
  push(@cmd_words, '--solid_nocall_strategy', 'LEAVE_READ_UNRECALIBRATED') if ($self->options('solid_nocall_strategy'));

  push(@cmd_words, '-R', $self->reference);
  push(@cmd_words, '-I', $self->input_bam);

  foreach my $vcf (@{$self->known_sites_files}) {
    push(@cmd_words, '-knownSites', $vcf);
  }
  push(@cmd_words, '-o', $covariates_file);
  my $cmd = join(' ', @cmd_words);

  $self->covariates_file($covariates_file);
  $self->execute_command_line ($cmd);

  return;
}


sub create_recalibrated_bam {
  my $self = shift;

  my $recalibrated_bam = $self->working_dir . '/'
        . $self->job_name. '.recal.bam';
  $recalibrated_bam =~ s{//}{/}g;

  my @cmd_words = ($self->java_exe, $self->jvm_args, '-jar');
  push(@cmd_words, $self->gatk_path . '/' . $self->jar_file);
  push(@cmd_words, '-T', 'PrintReads');
  push(@cmd_words, '-l', $self->options('log_level')) if ($self->options('log_level'));
  if ($self->options('intervals')) {
    foreach my $interval (split(/,/, $self->options('intervals'))) {
      push(@cmd_words, '-L', $interval);
    }
  }
  push(@cmd_words, '-allowPotentiallyMisencodedQuals') if ($self->options('allowPotentiallyMisencodedQuals'));
  push(@cmd_words, '--disable_bam_indexing');
  push(@cmd_words, '--solid_nocall_strategy', 'LEAVE_READ_UNRECALIBRATED') if ($self->options('solid_nocall_strategy'));
  push(@cmd_words, '-BQSR', $self->covariates_file);
  push(@cmd_words, '-R', $self->reference);
  push(@cmd_words, '-I', $self->input_bam);
  push(@cmd_words, '-o', $recalibrated_bam);

  my $cmd = join(' ', @cmd_words);

  $self->output_files($recalibrated_bam);
  $self->execute_command_line ($cmd);

  return;

}

sub covariates_file {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{covariates_file} = $arg;
    $self->created_files($arg);
  }
  return $self->{covariates_file};
}

1;
