
package ReseqTrack::HiveProcess::RunSamtools;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use ReseqTrack::Tools::RunSamtools;


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

    my @allowed_cmds = qw(merge sort index fix_and_calmd calmd fixmate sam_to_bam);
    throw("Don't recognise command $command. Acceptable commands are: @allowed_cmds")
      if (! grep {$command eq $_ } @allowed_cmds);

    my $samtools_object = ReseqTrack::Tools::RunSamtools->new(
      -input_files  => $bams,
      -program      => $program_file,
      -working_dir  => $output_dir,
      -job_name     => $job_name,
      -reference    => $self->param('reference'),
    );

    $samtools_object->run($command);
    $self->output_this_branch($command eq 'index' ? 'bai' : 'bam'  => $samtools_object->output_files);

}


1;

