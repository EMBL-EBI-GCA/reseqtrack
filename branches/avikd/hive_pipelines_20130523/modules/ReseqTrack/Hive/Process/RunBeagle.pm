package ReseqTrack::Hive::Process::RunBeagle;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::RunPhase::PhaseByBeagle;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);

=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    $self->param_required('vcf');

    my $vcfs = $self->file_param_to_flat_array('vcf');
 
    foreach my $vcf (@$vcfs) {
      check_file_exists($vcf);
    }
    
    my $fai = $self->param_required('fai');
    my $SQ_start = $self->param_required('SQ_start');
    my $SQ_end = $self->param_required('SQ_end');
    my $bp_start = $self->param_required('bp_start');
    my $bp_end = $self->param_required('bp_end');
    my $overlap = $self->param_is_defined('overlap') ? $self->param('overlap') : 0;
    
    my $slices = fai_to_slices(
          fai => $fai,
          SQ_start => $SQ_start, SQ_end => $SQ_end,
          bp_start => $bp_start, bp_end => $bp_end,
          );
    
    my $bed = $self->param_is_defined('bed') ? $self->param('bed') : undef;
         
    if (defined $bed) {
      $slices = bed_to_slices(bed => $bed, parent_slices => $slices);
    }
    
    throw("expected exactly one slice, not ". scalar @$slices) if scalar @$slices != 1;

    $slices->[0]->extend($overlap);
    

    my $run_phase = ReseqTrack::Tools::RunPhase::PhaseByBeagle->new(
          -input_files => $vcfs,
          -program => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
          -java_exe => $self->param_is_defined('java_exe') ? $self->param('java_exe') : undef,
          -jvm_args => $self->param_is_defined('jvm_args') ? $self->param('jvm_args') : undef,
          -working_dir => $self->output_dir,
          -reference_config => $self->param_is_defined('reference_config') ? $self->param('reference_config') : undef,
          -job_name => $self->job_name,
          -options => $self->param_is_defined('beagle_options') ? $self->param('beagle_options') : undef,
          -chrom => $slices->[0]->SQ_name ,
          -gl => $self->param_is_defined('gl') ? $self->param('gl') : 0,
          );
          
    if (!$slices->[0]->is_whole_SQ) {
      $run_phase->region_start($slices->[0]->start);
      $run_phase->region_end($slices->[0]->end);
    }

    $self->run_program($run_phase);

    $self->output_param('beagle_vcf', $run_phase->output_files->[0]);

}

1;

