
package ReseqTrack::Hive::Process::RunValidateBam;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunValidateBam;


sub param_defaults {
  return {
    program_file => undef,
  };
}

sub run {
    my ($self) = @_;

    $self->param_required('bam');
    my $bams = $self->param_as_array('bam');

    my $bam_validator = ReseqTrack::Tools::RunValidateBam->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -program      => $self->param('program_file'),
    );

    $self->run_program($bam_validator);

    my $output_files = ref($self->param('bam')) ? $bam_validator->output_files : $bam_validator->output_files->[0];
    $self->output_param('bas' => $output_files);

}


1;

