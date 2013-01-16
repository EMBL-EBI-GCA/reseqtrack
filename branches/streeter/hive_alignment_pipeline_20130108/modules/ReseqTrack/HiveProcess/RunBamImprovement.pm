
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
    my $type_output = $self->param('type_output') || die "'type_bam' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";
    my $java_exe = $self->param('java_exe');
    my $jvm_args = $self->param('jvm_args');
    my $gatk_dir = $self->param('gatk_dir');
    my $reference = $self->param('reference') || die "'reference' is an obligatory parameter";
    my $options = $self->param('gatk_module_options');
    my $branch_label = $self->param('branch_label') || die "'branch_label' is an obligatory parameter";
    my $command = $self->param('command') || die "'command' is an obligatory parameter";

    my $process_label = $self->param('process_label') || $command;

    throw ("Expecting one bam file") if ref($bams) eq 'ARRAY' && scalar @$bams != 1;

    my $bam = ref($bams) eq 'ARRAY' ? $bams->[0] : $bams;
    check_file_exists($bam);
    check_file_exists("$bam.bai");

    check_directory_exists($output_dir);
    my $output_bam = "$output_dir/$branch_label.$process_label.bam";

    system("touch $output_bam");

    $self->output_this_branch($type_output => $output_bam);
}


1;

