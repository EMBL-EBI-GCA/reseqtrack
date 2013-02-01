
package ReseqTrack::HiveProcess::BWA;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use ReseqTrack::Tools::RunAlignment::BWA;
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $fastqs = $self->param('fastq') || die "'fastq' is an obligatory parameter";
    my $run_id = $self->param('run_id') || die "'run_id' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";
    my $job_name = $self->param('job_name') or die "'job_name' is an obligatory parameter";

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $rmia = $db->get_RunMetaInfoAdaptor;
    my $rmi = $rmia->fetch_by_run_id($run_id);
    throw("did not get run_meta_info for $run_id") if !$rmi;

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

    $fastqs = ref($fastqs) eq 'ARRAY' ? $fastqs : [$fastqs];

    foreach my $fastq (@$fastqs) {
      check_file_exists($fastq);
    }

    $self->data_dbc->disconnect_when_inactive(1);

    my $run_alignment = ReseqTrack::Tools::RunAlignment::BWA->new(
          -input_files => $fastqs,
          -program => $self->param('program_file'),
          -samtools => $self->param('samtools'),
          -output_format => 'BAM',
          -working_dir => $output_dir,
          -reference => $self->param('reference'),
          -job_name => $job_name,
          -paired_length => $rmi->paired_length,
          -read_group_fields => \%read_group_fields,
          );
    $run_alignment->run;

    $self->output_this_branch('bam' => $run_alignment->output_files);
    $self->data_dbc->disconnect_when_inactive(0);

}

1;

