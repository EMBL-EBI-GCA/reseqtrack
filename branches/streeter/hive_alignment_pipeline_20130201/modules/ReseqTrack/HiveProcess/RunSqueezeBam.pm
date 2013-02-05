
package ReseqTrack::HiveProcess::RunSqueezeBam;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunBamSqueeze;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my ($self) = @_;

    my $bams = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $program_file = $self->param('program_file');

    $self->data_dbc->disconnect_when_inactive(1);

    my $bam_squeezer = ReseqTrack::Tools::RunBamSqueeze->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -program      => $program_file,
      -job_name     => $self->job_name,
      -rm_tag_types => $self->param('rm_tag_types'),
      -options      => {keepOQ => $self->param('rm_OQ_fields') ? 0 : 1,
                        keepDups => $self->param('rm_dups') ? 0 : 1},
    );
    $bam_squeezer->run;
    $self->output_this_branch('bam' => $bam_squeezer->output_files);
    $self->data_dbc->disconnect_when_inactive(0);
}


1;

