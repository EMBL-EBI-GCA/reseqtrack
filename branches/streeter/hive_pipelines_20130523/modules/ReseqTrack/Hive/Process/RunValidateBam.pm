
package ReseqTrack::Hive::Process::RunValidateBam;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunValidateBam;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my ($self) = @_;

    $self->param_required('bam');
    my $bams = $self->file_param_to_flat_array('bam');

    $self->data_dbc->disconnect_when_inactive(1);

    my $bam_validator = ReseqTrack::Tools::RunValidateBam->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -program      => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
    );

    $self->run_program($bam_validator);

    my $output_files = ref($self->param('bam')) ? $bam_validator->output_files : $bam_validator->output_files->[0];
    $self->output_param('bas' => $output_files);

}


1;

