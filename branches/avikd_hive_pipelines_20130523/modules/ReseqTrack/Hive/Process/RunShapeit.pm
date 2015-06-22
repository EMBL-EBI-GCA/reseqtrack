package ReseqTrack::Hive::Process::RunShapeit;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::RunPhase::PhaseByShapeit;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);
use ReseqTrack::Hive::Utils::ImputeSliceUtils qw( get_chrX );
=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    $self->param_required('gen');
    $self->param_required('samples');

    my $gens = $self->file_param_to_flat_array('gen');
    my $samples = $self->file_param_to_flat_array('samples');
    my $chrX_string = $self->param_is_defined('chrX_string') ? $self->param('chrX_string') : undef;
    
    throw("expecting single gen file") if(scalar @$gens > 1);
    throw("gen and samples file mismatch") if( scalar @$gens != scalar @$samples);
     
    foreach my $gen (@$gens) {
      check_file_exists($gen);
    }
    
    foreach my $sample (@$samples) {
      check_file_exists($sample);
    }
    
    my $chrX_PAR1_start = $self->param_is_defined('chrX_PAR1_start') ? $self->param('chrX_PAR1_start') : 1;
    my $chrX_PAR1_end = $self->param_is_defined('chrX_PAR1_end') ? $self->param('chrX_PAR1_end') : 2699500;
    my $chrX_PAR2_start = $self->param_is_defined('chrX_PAR2_start') ? $self->param('chrX_PAR2_start') : 154933000;
    my $chrX_PAR2_end = $self->param_is_defined('chrX_PAR2_end') ? $self->param('chrX_PAR2_end') : 155270000;
    
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
    
    my $chrX =  ( $slices->[0]->SQ_name eq $chrX_string ) ? 1 : undef;
    
    my $chrom = $slices->[0]->SQ_name;
    
    if($chrX){
      if($slices->[0]->is_whole_SQ){
        $chrom = "X";
      }
      else {
        $chrom = get_chrX(
               bp_start => $bp_start, 
               bp_end => $bp_end, 
               chrX_PAR1_start => $chrX_PAR1_start, 
               chrX_PAR1_end => $chrX_PAR1_end, 
               chrX_PAR2_start => $chrX_PAR2_start, 
               chrX_PAR2_end => $chrX_PAR2_end,
              );
        #if($slices->[0]->start <= $par_1_end && $slices->[0]->end <= $par_1_end ){
        #  $chrom = "X_PAR1";
        #}
        #elsif($slices->[0]->start > $par_1_end && $slices->[0]->end < $par_2_start){
        #  $chrom = "X_nonPAR";
        #}
        #elsif($slices->[0]->start => $par_2_start && $slices->[0]->end => $par_2_start){
        #  $chrom = "X_PAR2";
        #}
      }    
    }
    
    
    my $run_phase = ReseqTrack::Tools::RunPhase::PhaseByShapeit->new(
          -input_files => $gens,
          -input_samples => $samples, 
          -program => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
          -working_dir => $self->output_dir,
          -reference_config => $self->param_is_defined('reference_config') ? $self->param('reference_config') : undef,
          -job_name => $self->job_name,
          -options => $self->param_is_defined('shapeit_options') ? $self->param('shapeit_options') : undef,
          -chrom => $chrom,
          -phase_without_ref => $self->param_is_defined('phase_without_ref') ? $self->param('phase_without_ref') : undef,
          -exclude_strand_flip => $self->param_is_defined('exclude_strand_flip') ? $self->param('exclude_strand_flip') : undef,
          -no_check => $self->param_is_defined('no_check') ? $self->param('no_check') : undef,
          -chrX => $chrX,          
          );
          
    if (!$slices->[0]->is_whole_SQ) {
      $run_phase->region_start($slices->[0]->start);
      $run_phase->region_end($slices->[0]->end);
    }

    $self->run_program($run_phase);
    
    throw("expecting single haps")  if scalar @{$run_phase->output_haps_files} > 1;

    $self->output_param('shapeit_haps', $run_phase->output_haps_files);
    $self->output_param('shapeit_samples', $run_phase->output_samples_files);

}

1;    