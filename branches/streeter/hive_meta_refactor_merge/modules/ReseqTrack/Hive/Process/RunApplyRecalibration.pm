
package ReseqTrack::Hive::Process::RunApplyRecalibration;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::RunVariantCall::RunApplyRecalibration;


sub param_defaults {
  return {
    java_exe => undef,
    jvm_args => undef,
    gatk_dir => undef,
    options => undef,
  };
}

sub run {
    my ($self) = @_;

    $self->param_required('vcf');
    $self->param_required('recal_file');
    $self->param_required('tranches_file');
    my $reference = $self->param_required('reference');

    my $vcf_arr = $self->param_as_array('vcf');
    my $tranches_arr = $self->param_as_array('tranches_file');
    my $recal_arr = $self->param_as_array('recal_file');

    throw("Expecting one vcf file") if scalar @$vcf_arr != 1;
    throw("Expecting one tranches file") if scalar @$tranches_arr != 1;
    throw("Expecting one recal file") if scalar @$recal_arr != 1;

    my $gatk_object = ReseqTrack::Tools::RunVariantCall::RunApplyRecalibration->new(
      -input_vcf    => $vcf_arr->[0],
      -tranches_file => $tranches_arr->[0],
      -recal_file   => $recal_arr->[0],
      -working_dir  => $self->output_dir,
      -reference    => $reference,
      -job_name     => $self->job_name,
      -java_exe     => $self->param('java_exe'),
      -jvm_args     => $self->param('jvm_args'),
      -gatk_path    => $self->param('gatk_dir'),
      -options      => $self->param('options'),
    );

    $self->run_program($gatk_object);

    $self->output_param('vcf', $gatk_object->output_files->[0]);
}


1;

