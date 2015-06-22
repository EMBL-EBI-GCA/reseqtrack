
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
    my $program_file = $self->param('program_file');
    my $options = $self->param('samtools_options') || {};

    my @allowed_cmds = qw(merge sort index fix_and_calmd calmd fixmate sam_to_bam);
    throw("Don't recognise command $command. Acceptable commands are: @allowed_cmds")
      if (! grep {$command eq $_ } @allowed_cmds);

    my @allowed_options = keys %{&ReseqTrack::Tools::RunSamtools::DEFAULT_OPTIONS};
    foreach my $option (keys %$options) {
      throw("Don't recognise option $option. Acceptable options are: @allowed_options")
        if (! grep {$option eq $_ } @allowed_options);
    }


    my $samtools_object = ReseqTrack::Tools::RunSamtools->new(
      -input_files  => $bams,
      -program      => $program_file,
      -working_dir  => $self->output_dir,
      -job_name     => $self->job_name,
      -reference    => $self->param('reference'),
      -options      => $options,
    );

    $self->data_dbc->disconnect_when_inactive(1);
    eval{$samtools_object->run($command);};
    my $msg_thrown = $@;
    $self->data_dbc->disconnect_when_inactive(0);
    die $msg_thrown if $msg_thrown;

    $self->output_this_branch($command eq 'index' ? 'bai' : 'bam'  => $samtools_object->output_files);

}


1;

