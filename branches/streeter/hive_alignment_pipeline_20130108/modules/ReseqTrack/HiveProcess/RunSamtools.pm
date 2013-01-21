
package ReseqTrack::HiveProcess::RunSamtools;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $bams = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $command = $self->param('command') || die "'command' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') or $command eq 'index' or die "'output_dir' is an obligatory parameter";
    my $job_name = $self->param('job_name') or die "'job_name' is an obligatory parameter";
    my $program_file = $self->param('program_file');

    $bams = ref($bams) eq 'ARRAY' ? $bams : [$bams];

    foreach my $bam (@$bams) {
      check_file_exists($bam);

      if ($command eq 'index') {
        my $bai = "$bam.bai";
        system("touch $bai");
        $self->output_this_branch('bai' => $bai);
      }
      else {
        check_directory_exists($output_dir);
        my $output_bam = "$output_dir/$job_name.bam";

        system("touch $output_bam");
        $self->output_this_branch('bam' => $output_bam);
      }
    }

}


1;

