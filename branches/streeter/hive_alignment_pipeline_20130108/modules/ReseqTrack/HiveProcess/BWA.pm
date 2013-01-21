
package ReseqTrack::HiveProcess::BWA;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
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
    my $program_file = $self->param('program_file');

    $fastqs = ref($fastqs) eq 'ARRAY' ? $fastqs : [$fastqs];

    foreach my $fastq (@$fastqs) {
      check_file_exists($fastq);
    }

    check_directory_exists($output_dir);

    my $bam = "$output_dir/$job_name.bam";

    system("touch $bam");

    $self->output_this_branch('bam' => $bam);

}

1;

