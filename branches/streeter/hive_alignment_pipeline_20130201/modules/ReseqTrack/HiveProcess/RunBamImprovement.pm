
package ReseqTrack::HiveProcess::RunBamImprovement;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my ($self) = @_;

    my $bams = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $java_exe = $self->param('java_exe');
    my $jvm_args = $self->param('jvm_args');
    my $gatk_dir = $self->param('gatk_dir');
    my $reference = $self->param('reference') || die "'reference' is an obligatory parameter";
    my $options = $self->param('gatk_module_options');
    my $command = $self->param('command') || die "'command' is an obligatory parameter";

    throw("Expecting one bam file") if ref($bams) eq 'ARRAY' && scalar @$bams != 1;
    throw("Expecting either 'realign' or 'recalibrate' for command")
        if ($command ne 'realign' && $command ne 'recalibrate');
    my $module_name = $command eq 'realign' ? 'IndelRealigner' : 'QualityScoreRecalibrator';
    eval{require "ReseqTrack/Tools/GATKTools/$module_name.pm";};
    throw("cannot load $module_name: $@") if $@;
    my $gatk_module = "ReseqTrack::Tools::GATKTools::$module_name";

    $self->data_dbc->disconnect_when_inactive(1);

    my $gatk_object = $gatk_module->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -reference    => $reference,
      -job_name     => $self->job_name,
      -java_exe     => $java_exe,
      -jvm_args     => $jvm_args,
      -gatk_path    => $gatk_dir,
      -options      => $options,
      -known_sites_files => $self->param('known_sites_vcf'),
    );
    $gatk_object->run($command);
    $self->output_this_branch('bam' => $gatk_object->output_bam);
    $self->data_dbc->disconnect_when_inactive(0);
}


1;

