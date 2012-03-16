=pod

=head1 NAME

ReseqTrack::Tools::GATKTools::QualityScoreRecalibrator

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

package ReseqTrack::Tools::GATKTools::QualityScoreRecalibrator;

use strict;
use warnings;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::GATKTools;

use vars qw(@ISA);
use Data::Dumper;

@ISA = qw(ReseqTrack::Tools::GATKTools);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);


        #setting defaults
        if (!$self->jar_file) {
          $self->jar_file ("GenomeAnalysisTK.jar");
        }
        if (!$self->options('CountCovariates')) {
          my $cc_options = "-l INFO -L '1;2;3;4;5;6;7;8;9;10;11;12;";
          $cc_options .= "13;14;15;16;17;18;19;20;21;22;X;Y;MT' ";
          $cc_options .= "-cov ReadGroupCovariate -cov QualityScoreCovariate ";
          $cc_options .= "-cov CycleCovariate -cov DinucCovariate";
          $self->options ('CountCovariates', $cc_options);
        }

        if (!$self->options('TableRecalibration')) {
          my $tr_options = "-l INFO -compress 0 --disable_bam_indexing";
          $self->options ('TableRecalibration', $tr_options);
        }

	return $self;
}


sub run_program {
	my $self = shift;

        throw "no input bam" if (!$self->input_bam);
        throw "\nNo CountCovariate options\n"
            if ( ! $self->{options}->{'CountCovariates'});
        throw "\nNo TableRecalibration options\n"
            if ( ! $self->{options}->{'TableRecalibration'});
        warn "No known sites files"
            if (! @{$self->known_sites_files});
        check_file_exists($_) foreach (@{$self->known_sites_files});
        $self->check_jar_file_exists;
        check_file_exists($self->reference);


	$self->create_recalibration_table();
	$self->create_recalibrated_bam();

	return;
}


sub create_recalibration_table {
  my $self = shift;

  my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
  $cmd .= $self->gatk_path ;
  $cmd .= "\/";
  $cmd .= $self->jar_file;
  $cmd .= " -T CountCovariates ";

  $cmd .= $self->{options}->{'CountCovariates'};

  $cmd .= " -R " . $self->reference;
  $cmd .= " -I " . $self->input_bam;

  foreach my $vcf (@{$self->known_sites_files}) {
    $cmd .= " -knownSites $vcf";
  }

  my $covariates_file = $self->working_dir . '/'
        . $self->job_name . '.bam.recal_data.csv';
  $covariates_file =~ s{//}{/}g;

  $cmd .= " -recalFile $covariates_file";

  $self->covariates_file($covariates_file);
  $self->execute_command_line ($cmd);


  return;
}


sub create_recalibrated_bam {
  my $self = shift;

  my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
  $cmd .= $self->gatk_path . "\/" . $self->jar_file;
  $cmd .= " -T TableRecalibration ";
  $cmd .= $self->{options}->{'TableRecalibration'};

  $cmd .= " -recalFile " . $self->covariates_file . " ";
  $cmd .= " -R " . $self->reference . " ";
  $cmd .= " -I " . $self->input_bam;

  my $recalibrated_bam = $self->working_dir . '/'
        . $self->job_name. '.recal.bam';
  $cmd .= " -o " . $recalibrated_bam;

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
