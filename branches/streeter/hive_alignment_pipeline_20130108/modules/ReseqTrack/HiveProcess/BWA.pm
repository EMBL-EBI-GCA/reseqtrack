
package ReseqTrack::HiveProcess::BWA;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists);
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $fastq = $self->param('fastq') || die "'fastq' is an obligatory parameter";
    my $run_id = $self->param('run_id') || die "'run_id' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";
    my $branch_label = $self->param('branch_label') || die "'branch_label' is an obligatory parameter";
    my $type_bam = $self->param('type_bam') || die "'type_bam' is an obligatory parameter";
    my $program_file = $self->param('program_file');

    my $process_label = $self->param('process_label') || 'bwa';

    check_directory_exists($output_dir);

    my $bam = "$output_dir/$branch_label.$process_label.bam";

    system("touch $bam");

    $self->output_this_branch($type_bam => $bam);

}

1;

