
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

    my $run = $db->get_RunAdaptor->fetch_by_dbID($run_id);
    throw('did not get run with id '.$run_id) if !$run;
    my $sample = $db->get_SampleAdaptor->fetch_by_dbID($run->sample_id);
    throw('did not get sample with id '.$run->sample_id) if !$sample;
    my $experiment = $db->get_ExperimentAdaptor->fetch_by_dbID($run->experiment_id);
    throw('did not get experiment with id '.$run->experiment_id) if !$experiment;
    my $study = $db->get_StudyAdaptor($experiment->study_id)->fetch_by_dbID($experiment->study_id);
    throw('did not get study with id '.$experiment->study_id) if !$study;

    $db->dbc->disconnect_when_inactive(1);

    my %read_group_fields = (
      ID => $run->source_id,
      CN => $run->center_name,
      LB => $experiment->library_name,
      PI => $experiment->paired_nominal_length,
      SM => $sample->sample_alias,
      DS => $study->source_id,
      PU => $run->run_alias,
      PL => $experiment->instrument_platform,
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
          -paired_length => $experiment->paired_nominal_length,
          -read_group_fields => \%read_group_fields,
          );

    $self->run_program($run_alignment);

    $self->output_param('bam', $run_alignment->output_files->[0]);

}

1;

