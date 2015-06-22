package ReseqTrack::Hive::Process::RunVcfToImpute;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunVcfToImpute;
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);

=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    $self->param_required('vcf');
    my $vcfs = $self->file_param_to_flat_array('vcf');
    
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


    my $vcfToImpute_object = ReseqTrack::Tools::RunVcfToImpute->new(
      -input_files  => $vcfs,
      -working_dir  => $self->output_dir,
      -job_name     => $self->job_name,
      -program      => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
      -options => $self->param_is_defined('options') ? $self->param('options') : undef,
      -chrom => $slices->[0]->SQ_name,
      -region_start => $bp_start,
      -region_end => $bp_end,
    );

    $self->run_program($vcfToImpute_object);

    my $output_gen = $vcfToImpute_object->output_gen_files;
    my $output_samples = $vcfToImpute_object->output_samples_files;
    
    if (@$output_gen ==1) {
      $output_gen = $output_gen->[0];
    }
    if (@$output_samples ==1) {
      $output_samples = $output_samples->[0];
    }

    $self->output_param('gen', $output_gen);
    $self->output_param('samples', $output_samples);

}


1;