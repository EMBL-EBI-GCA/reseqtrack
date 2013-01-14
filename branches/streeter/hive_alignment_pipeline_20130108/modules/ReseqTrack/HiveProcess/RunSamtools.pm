
package ReseqTrack::HiveProcess::RunSamtools;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $bam = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $type_bam = $self->param('type_bam') || die "'type_bam' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";
    my $program_file = $self->param('program_file');
    my $branch_label = $self->param('branch_label') || die "'branch_label' is an obligatory parameter";
    my $command = $self->param('command') || die "'command' is an obligatory parameter";

    my $process_label = $self->param('process_label') || $command;

    check_directory_exists($output_dir);
    my $bam = "$output_dir/$branch_label.$process_label.bam";

    system("touch $bam");

    $self->output_this_branch($type_bam => $bam);
}


1;

