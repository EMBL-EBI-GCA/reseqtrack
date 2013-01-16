
package ReseqTrack::HiveProcess::RunPicard;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $bams = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $command = $self->param('command') || die "'command' is an obligatory parameter";
    my $type_output = $self->param('type_output') or die "'type_output' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') or die "'output_dir' is an obligatory parameter";
    my $branch_label = $self->param('branch_label') or die "'branch_label' is an obligatory parameter";
    my $java_exe = $self->param('java_exe');
    my $jvm_args = $self->param('jvm_args');
    my $picard_dir = $self->param('picard_dir');

    $bams = ref($bams) eq 'ARRAY' ? $bams : [$bams];

    foreach my $bam (@$bams) {
      check_file_exists($bam);

      my $process_label = $self->param('process_label') || $command;

      check_directory_exists($output_dir);
      my $bam = "$output_dir/$branch_label.$process_label.bam";

      system("touch $bam");
      $self->output_this_branch($type_output => $bam);
    }

}


1;

