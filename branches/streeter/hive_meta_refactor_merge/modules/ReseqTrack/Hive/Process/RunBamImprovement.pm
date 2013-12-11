
package ReseqTrack::Hive::Process::RunBamImprovement;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);


sub param_defaults {
  return {
    java_exe => undef,
    jvm_args => undef,
    gatk_dir => undef,
    gatk_module_options => undef,
    known_sites_vcf => undef,
    intervals_file => undef,
  };
}


sub run {
    my ($self) = @_;

    $self->param_required('bam');
    my $bams = $self->param_as_array('bam');
    my $reference = $self->param_required('reference');
    my $command = $self->param_required('command');

    throw("Expecting one bam file") if scalar @$bams != 1;
    throw("Expecting either 'realign' or 'recalibrate' for command")
        if ($command ne 'realign' && $command ne 'recalibrate');
    my $module_name = $command eq 'realign' ? 'IndelRealigner' : 'QualityScoreRecalibrator';
    eval{require "ReseqTrack/Tools/GATKTools/$module_name.pm";};
    throw("cannot load $module_name: $@") if $@;
    my $gatk_module = "ReseqTrack::Tools::GATKTools::$module_name";


    my $gatk_object = $gatk_module->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -reference    => $reference,
      -job_name     => $self->job_name,
      -java_exe     => $self->param('java_exe'),
      -jvm_args     => $self->param('jvm_args'),
      -gatk_path    => $self->param('gatk_dir'),
      -options      => $self->param('gatk_module_options'),
      -known_sites_files => $self->param('known_sites_vcf'),
    );
    if ($command eq 'realign' && defined $self->param('intervals_file')) {
      $gatk_object->intervals_file($self->param('intervals_file'));
    }

    $self->run_program($gatk_object, $command);

    $self->output_param('bam', $gatk_object->output_bam);
}


1;

