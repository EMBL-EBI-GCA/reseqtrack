
package ReseqTrack::Hive::Process::RunVariantRecalibrator;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my ($self) = @_;

    $self->param_required('vcf');
    $self->param_required('resources');
    my $vcfs = $self->file_param_to_flat_array('vcf');
    my $reference = $self->param_required('reference');

    throw("Expecting one vcf file") if scalar @$vcfs != 1;


    my %resources;
    while( my ($name, $string) = each %{$self->param('resources')}) {
      my ($resource_options, $resource_file) = split(/\s+/, $string);
      throw("resource string should be e.g. 'known=false,training=true,truth=false,prior=12.0 /path/to/file': $string")
          if ! $resource_file;
      $resources{$name.','.$resource_options} = $resource_file;
    }

    my $gatk_object = ReseqTrack::Tools::RunVariantCall::RunVariantRecalibrator->new(
      -input_files  => $vcfs,
      -working_dir  => $self->output_dir,
      -reference    => $reference,
      -job_name     => $self->job_name,
      -java_exe     => $self->param_is_defined('java_exe') ? $self->param('java_exe') : undef,
      -jvm_args     => $self->param_is_defined('jvm_args') ? $self->param('jvm_args') : undef,
      -gatk_path    => $self->param_is_defined('gatk_dir') ? $self->param('gatk_dir') : undef,
      -options      => $self->param_is_defined('options') ? $self->param('options') : undef,
      -annotations  => $self->param_is_defined('annotations') ? $self->param('annotations') : undef,
      -resources    => \%resources
    );

    $self->run_program($gatk_object);

    $self->output_param('recal_file', $gatk_object->output_recal_file);
    $self->output_param('tranches_file', $gatk_object->output_tranches_file);
    $self->output_param('tranches_pdf', $gatk_object->output_pdf_file);
}


1;

