
package ReseqTrack::Hive::Process::RunPicard;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunPicard;
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    $self->param_required('bam');
    my $bams = $self->get_param_values('bam');
    my $command = $self->param_required('command');

    my @allowed_cmds = ReseqTrack::Tools::RunPicard->get_valid_commands;
    throw( "Don't recognise command $command. Acceptable commands are: @allowed_cmds")
      if ( !grep { $command eq $_ } @allowed_cmds );


    my $picard_object = ReseqTrack::Tools::RunPicard->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -job_name     => $self->job_name,
      -java_exe     => $self->param_is_defined('java_exe') ? $self->param('java_exe') : undef,
      -jvm_options  => $self->param_is_defined('jvm_args') ? $self->param('jvm_args') : undef,
      -picard_dir     => $self->param_is_defined('picard_dir') ? $self->param('picard_dir') : undef,
      -options      => {validation_stringency => 'SILENT'},
      -create_index => $self->param_is_defined('create_index') ? $self->param('create_index') : undef,
      -keep_metrics => 0,
    );

    $self->run_program($picard_object, $command);

    my $output_bams = $picard_object->output_bam_files;
    my $output_bais = $picard_object->output_bai_files;
    if (@$output_bams ==1) {
      $output_bams = $output_bams->[0];
    }
    if (@$output_bais ==1) {
      $output_bais = $output_bais->[0];
    }

    $self->output_param('bam', $output_bams);
    $self->output_param('bai', $output_bais);

}


1;

