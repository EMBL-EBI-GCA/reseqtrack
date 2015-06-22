
package ReseqTrack::Hive::Process::RunVcfTools;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use ReseqTrack::Tools::RunVcfTools;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    $self->param_required('vcf');
    my $vcfs = $self->file_param_to_flat_array('vcf');
    my $command = $self->param_required('command');
    my $vcftools_dir = $self->param_required('vcftools_dir');
    my $reference_index = $self->param_required('reference_index');

    my @allowed_cmds = qw( merge concat sort vcftools );
    throw("Don't recognise command $command. Acceptable commands are: @allowed_cmds")
      if (! grep {$command eq $_ } @allowed_cmds);

#     my $options = $self->param_is_defined('vcftools_options') ? $self->param('vcftools_options') : undef;
#     my @allowed_options = keys %{&ReseqTrack::Tools::RunVcfTools::DEFAULT_OPTIONS};
#     foreach my $option (keys %$options) {
#       throw("Don't recognise option $option. Acceptable options are: @allowed_options")
#         if (! grep {$option eq $_ } @allowed_options);
#     }


    my $vcftools_object = ReseqTrack::Tools::RunVcfTools->new(
      -input_files  => $vcfs,
      -program      => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
      -working_dir  => $self->output_dir,
      -job_name     => $self->job_name,
      -vcftools_dir    => $vcftools_dir,
      -tabix           => $self->param_is_defined('tabix') ? $self->param('tabix') : undef,
      -bgzip           => $self->param_is_defined('bgzip') ? $self->param('bgzip') : undef,
      -create_index    => $self->param_is_defined('create_index') && $self->param('tabix') ? 1 : 0,
      -reference_index => $reference_index,
    );

    $self->run_program($vcftools_object, $command);

    my $output_vcfs = $vcftools_object->output_vcf_files;
    $self->output_param('vcf',$output_vcfs);
    
    if($self->param_is_defined('create_index')) {
      my $output_tbi = $vcftools_object->output_tbi_files;
      $self->output_param('tbi',$output_vcfs);
    }
    
    

    

}


1;