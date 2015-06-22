package ReseqTrack::Hive::Process::RunImpute2;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::RunImpute::ImputeByImpute2;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);
use ReseqTrack::Hive::Utils::ImputeSliceUtils qw( get_chrX );

=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    $self->param_required('haps');
    $self->param_required('samples');
    $self->param_required('reference_config');

    my $haps = $self->file_param_to_flat_array('haps');
    my $samples = $self->file_param_to_flat_array('samples');
    
    my $chrX_string = $self->param_is_defined('chrX_string') ? $self->param('chrX_string') : undef;
    my $use_phased = $self->param_is_defined('use_phased') ? $self->param('use_phased') : undef;
    my $phase_result = $self->param_is_defined('phase_result') ? $self->param('phase_result') : undef;
    my $reference_config = $self->param('reference_config');
    
    
    
    throw("expecting single haps file") if(scalar @$haps > 1);
    throw("haps and samples file mismatch") if( scalar @$haps != scalar @$samples);
    throw("phased input required for imputation") if($use_phased != 1);
    
    foreach my $hap (@$haps) {
      
      check_file_exists($hap);
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
        #  $chrom = "non_PAR";
        #}
        #elsif($slices->[0]->start => $par_2_start && $slices->[0]->end => $par_2_start){
        #  $chrom = "X_PAR2";
        #}
      }    
    }
    
    my $run_impute = ReseqTrack::Tools::RunImpute::ImputeByImpute2->new(
          -input_files => $haps,
          -input_samples => $samples, 
          -program => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
          -working_dir => $self->output_dir,
          -reference_config => $reference_config,
          -job_name => $self->job_name,
          -options => $self->param_is_defined('impute2_options') ? $self->param('impute2_options') : undef,
          -chrom => $chrom,
          -use_phased => $use_phased,
          -phase_result => $phase_result,
          -chrX => $chrX, 
          -region_start => $bp_start,
          -region_end => $bp_end,
          -keep_impute => $self->param_is_defined('keep_impute') ? $self->param('keep_impute') : undef,
          -keep_samples => $self->param_is_defined('keep_samples') ? $self->param('keep_samples') : undef,
          -keep_haps => $self->param_is_defined('keep_haps') ? $self->param('keep_haps') : undef,
          -keep_info => $self->param_is_defined('keep_info') ? $self->param('keep_info') : undef,
          -keep_allele_probs => $self->param_is_defined('keep_allele_probs') ? $self->param('keep_allele_probs') : undef,
          -keep_diplotype_ordering => $self->param_is_defined('keep_diplotype_ordering') ? $self->param('keep_diplotype_ordering') : undef,
          -keep_info_by_sample => $self->param_is_defined('keep_info_by_sample') ? $self->param('keep_info_by_sample') : undef,
          -keep_warnings => $self->param_is_defined('keep_warnings') ? $self->param('keep_warnings') : undef,
          -keep_summary => $self->param_is_defined('keep_summary') ? $self->param('keep_summary') : undef,
          );
          
    

    $self->run_program($run_impute);


    $self->output_param('impute',$run_impute->output_impute_files) if $self->param('keep_impute');
    $self->output_param('impute_haps',$run_impute->output_haps_files) if $phase_result && $self->param('keep_haps');
    $self->output_param('impute_allele_probs',$run_impute->output_allele_probs_files) if $self->param('keep_allele_probs') ;
    $self->output_param('impute_info',$run_impute->output_info_files) if $self->param('keep_info');
    $self->output_param('impute_diplotype_ordering',$run_impute->output_diplotype_ordering_files) if $self->param('keep_diplotype_ordering');
    $self->output_param('impute_info_by_sample',$run_impute->output_info_by_sample_files) if $self->param('keep_info_by_sample');
    $self->output_param('impute_samples',$run_impute->output_samples_files) if $self->param('keep_samples');
    $self->output_param('impute_summary',$run_impute->output_summary_files) if $self->param('keep_summary');
    $self->output_param('impute_warnings',$run_impute->output_warnings_files) if $self->param('keep_warnings');

}

1;    