
package ReseqTrack::HiveProcess::RunSqueezeBam;

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
    my $branch_label = $self->param('branch_label') || die "'branch_label' is an obligatory parameter";
    my $program_file = $self->param('program_file');
    my $rm_OQ_files = $self->param('rm_OQ_fields');
    my $rm_dups = $self->param('rm_dups');
    my $rm_tag_type = $self->param('rm_tag_type');

    my $process_label = $self->param('process_label') || 'squeezed';

    throw ("Expecting one bam file") if ref($bams) eq 'ARRAY' && scalar @$bams != 1;

    my $bam = ref($bams) eq 'ARRAY' ? $bams->[0] : $bams;
    check_file_exists($bam);

    check_directory_exists($output_dir);
    my $output_bam = "$output_dir/$branch_label.$process_label.bam";

    system("touch $output_bam");

    $self->output_this_branch($type_output => $output_bam);
}


1;

