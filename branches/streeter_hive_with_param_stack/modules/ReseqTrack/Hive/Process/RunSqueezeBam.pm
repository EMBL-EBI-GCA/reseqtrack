
package ReseqTrack::Hive::Process::RunSqueezeBam;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunBamSqueeze;

sub param_defaults {
  return {
    program_file => undef,
    rm_tag_types => undef,
    rm_OQ_fields => 1,
    rum_dups => 1,
  };
}


sub run {
    my ($self) = @_;

    $self->param_required('bam');
    my $bams = $self->file_param_to_flat_array('bam');
    throw('too many bam files: '. join(' ', @$bams)) if @$bams !=1;

    my $bam_squeezer = ReseqTrack::Tools::RunBamSqueeze->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -program      => $self->param('program_file'),
      -job_name     => $self->job_name,
      -rm_tag_types => $self->param('rm_tag_types'),
      -options      => {keepOQ => $self->param('rm_OQ_fields') ? 0 : 1,
                        keepDups => $self->param('rm_dups') ? 0 : 1},
    );

    $self->run_program($bam_squeezer);
    $self->output_param('bam' => $bam_squeezer->output_files->[0]);
}


1;

