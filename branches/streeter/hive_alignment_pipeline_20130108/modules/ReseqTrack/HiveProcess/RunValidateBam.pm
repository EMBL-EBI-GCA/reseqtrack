
package ReseqTrack::HiveProcess::RunValidateBam;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunValidateBam;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my ($self) = @_;

    my $bams = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";
    my $program_file = $self->param('program_file');

    $self->data_dbc->disconnect_when_inactive(1);

    my $bam_validator = ReseqTrack::Tools::RunValidateBam->new(
      -input_files  => $bams,
      -working_dir  => $output_dir,
      -program      => $program_file,
    );
    $bam_validator->run;
    $self->output_this_branch('bas' => $bam_validator->output_files);
    $self->data_dbc->disconnect_when_inactive(0);
}


1;

