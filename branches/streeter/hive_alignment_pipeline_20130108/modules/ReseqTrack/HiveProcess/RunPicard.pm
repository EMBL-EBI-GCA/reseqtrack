
package ReseqTrack::HiveProcess::RunPicard;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunPicard;
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $bams = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $command = $self->param('command') || die "'command' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') or die "'output_dir' is an obligatory parameter";
    my $job_name = $self->param('job_name') or die "'job_name' is an obligatory parameter";
    my $java_exe = $self->param('java_exe');
    my $jvm_args = $self->param('jvm_args');
    my $picard_dir = $self->param('picard_dir');

    my @allowed_cmds = ReseqTrack::Tools::RunPicard->get_valid_commands;
    throw( "Don't recognise command $command. Acceptable commands are: @allowed_cmds")
      if ( !grep { $command eq $_ } @allowed_cmds );

    my $picard_object = ReseqTrack::Tools::RunPicard->new(
      -input_files  => $bams,
      -working_dir  => $output_dir,
      -job_name     => $job_name,
      -java_exe     => $java_exe,
      -jvm_options  => $jvm_args,
      -picard_dir   => $picard_dir,
      -options      => {validation_stringency => 'SILENT'},
      -create_index => $self->param('create_index'),
      -keep_metrics => 0,
    );

    $picard_object->run($command);
    $self->output_this_branch('bam' => $picard_object->output_bam_files,
                              'bai' => $picard_object->output_bai_files);

}


1;

