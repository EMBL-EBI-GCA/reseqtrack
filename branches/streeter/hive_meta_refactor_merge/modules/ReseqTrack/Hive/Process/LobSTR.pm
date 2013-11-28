
package ReseqTrack::Hive::Process::LobSTR;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::RunAlignment::LobSTR;
use ReseqTrack::Tools::Exception qw(throw);

sub param_defaults {
  return {
    program_file => undef,
    options => {},
    library_name => undef,
    sample_source_id => undef,
  };
}

sub run {
    my $self = shift @_;

    $self->param_required('fastq');
    my $ref_index_prefix = $self->param_required('ref_index_prefix');

    my $fastqs = $self->file_param_to_flat_array('fastq');

    my %read_group_fields = (
      LB => $self->param('library_name'),
      SM => $self->param('sample_source_id'),
    );

    foreach my $fastq (@$fastqs) {
      check_file_exists($fastq);
    }


    my $run_alignment = ReseqTrack::Tools::RunAlignment::LobSTR->new(
          -input_files => $fastqs,
          -program => $self->param('program_file'),
          -working_dir => $self->output_dir,
          -job_name => $self->job_name,
          -ref_index_prefix => $ref_index_prefix,
          -read_group_fields => \%read_group_fields,
          -options => $self->param('options'),
          );

    $self->run_program($run_alignment);

    $self->output_param('bam', $run_alignment->output_files->[0]);

}

1;

