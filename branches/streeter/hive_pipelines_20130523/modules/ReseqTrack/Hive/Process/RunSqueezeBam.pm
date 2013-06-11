
package ReseqTrack::Hive::Process::RunSqueezeBam;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunBamSqueeze;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my ($self) = @_;

    $self->param_required('bam');
    my $bams = $self->get_param_values('bam');
    throw('too many bam files: '. join(' ', @$bams)) if @$bams !=1;

    my $bam_squeezer = ReseqTrack::Tools::RunBamSqueeze->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -program      => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
      -job_name     => $self->job_name,
      -rm_tag_types => $self->param_is_defined('rm_tag_types') ? $self->param('rm_tag_types') : undef,
      -options      => {keepOQ => $self->param_is_defined('rm_OQ_fields') && $self->param('rm_OQ_fields') ? 0 : 1,
                        keepDups => $self->param_is_defined('rm_dups') && $self->param('rm_dups') ? 0 : 1},
    );

    $self->run_program($bam_squeezer);
    $self->output_param('bam' => $bam_squeezer->output_files->[0]);
}


1;

