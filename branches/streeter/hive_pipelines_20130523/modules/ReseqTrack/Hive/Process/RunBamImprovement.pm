
package ReseqTrack::Hive::Process::RunBamImprovement;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my ($self) = @_;

    $self->param_required('bam');
    my $bams = $self->get_param_values('bam');
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
      -java_exe     => $self->param_is_defined('java_exe') ? $self->param('java_exe') : undef,
      -jvm_args     => $self->param_is_defined('jvm_args') ? $self->param('jvm_args') : undef,
      -gatk_path    => $self->param_is_defined('gatk_dir') ? $self->param('gatk_dir') : undef,
      -options      => $self->param_is_defined('gatk_module_options') ? $self->param('gatk_module_options') : undef,
      -known_sites_files => $self->param_is_defined('known_sites_vcf') ? $self->param('known_sites_vcf') : undef,
    );
    if ($command eq 'realign' && $self->param_is_defined('intervals_file')) {
      $gatk_object->intervals_file($self->param('intervals_file'));
    }

    $self->run_program($gatk_object, $command);

    $self->output_param('bam', $gatk_object->output_bam);
}


1;

