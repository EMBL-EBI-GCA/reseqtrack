
package ReseqTrack::Hive::Process::RunSamtools;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use ReseqTrack::Tools::RunSamtools;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    $self->param_required('bam');
    my $bams = $self->file_param_to_flat_array('bam');
    my $command = $self->param_required('command');

    my @allowed_cmds = qw(merge sort index fix_and_calmd calmd fixmate sam_to_bam);
    throw("Don't recognise command $command. Acceptable commands are: @allowed_cmds")
      if (! grep {$command eq $_ } @allowed_cmds);

    my $options = $self->param_is_defined('samtools_options') ? $self->param('samtools_options') : {};
    my @allowed_options = keys %{&ReseqTrack::Tools::RunSamtools::DEFAULT_OPTIONS};
    foreach my $option (keys %$options) {
      throw("Don't recognise option $option. Acceptable options are: @allowed_options")
        if (! grep {$option eq $_ } @allowed_options);
    }


    my $samtools_object = ReseqTrack::Tools::RunSamtools->new(
      -input_files  => $bams,
      -program      => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
      -working_dir  => $self->output_dir,
      -job_name     => $self->job_name,
      -reference    => $self->param_is_defined('reference') ? $self->param('reference') : undef,
      -options      => $options,
    );

    $self->run_program($samtools_object, $command);

    my $output_files = $samtools_object->output_files;
    if (!ref($self->param('bam')) || $command eq 'merge') {
      $output_files = $output_files->[0];
    }

    $self->output_param($command eq 'index' ? 'bai' : 'bam'  => $output_files);

}


1;

