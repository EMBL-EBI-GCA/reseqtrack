
package ReseqTrack::Hive::Process::BWA;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::RunAlignment::BWA;
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    $self->param_required('fastq');
    my $run_id = $self->param_required('run_id');
    my $reference = $self->param_required('reference');

    my $fastqs = $self->file_param_to_flat_array('fastq');


    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $rmia = $db->get_RunMetaInfoAdaptor;
    my $rmi = $rmia->fetch_by_run_id($run_id);
    throw("did not get run_meta_info for $run_id") if !$rmi;
    $db->dbc->disconnect_when_inactive(1);

    my %read_group_fields = (
      ID => $rmi->run_id,
      CN => $rmi->center_name,
      LB => $rmi->library_name,
      PI => $rmi->paired_length,
      SM => $rmi->sample_name,
      DS => $rmi->study_id,
      PU => $rmi->run_name,
      PL => $rmi->instrument_platform,
    );

    foreach my $fastq (@$fastqs) {
      check_file_exists($fastq);
    }


    my $run_alignment = ReseqTrack::Tools::RunAlignment::BWA->new(
          -input_files => $fastqs,
          -program => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
          -samtools => $self->param_is_defined('samtools') ? $self->param('samtools') : undef,
          -output_format => 'BAM',
          -working_dir => $self->output_dir,
          -reference => $reference,
          -job_name => $self->job_name,
          -paired_length => $rmi->paired_length,
          -read_group_fields => \%read_group_fields,
          );

    $self->run_program($run_alignment);

    $self->output_param('bam', $run_alignment->output_files->[0]);

}

1;

