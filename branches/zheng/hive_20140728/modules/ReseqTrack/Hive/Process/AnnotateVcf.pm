package ReseqTrack::Hive::Process::AnnotateVcf;

use strict;
use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::AddReadDpAddGlCalculateAF;
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw(fai_to_slices);

sub param_defaults {
  return {
	overlap => 0,
    save_files_from_deletion => 1,  ### FIXME, check this to 0
  };
}

sub run {
    my $self = shift @_;

    my $fai = $self->param_required('fai');
    my $SQ_start = $self->param_required('SQ_start');
    my $SQ_end = $self->param_required('SQ_end');
    my $bp_start = $self->param_required('bp_start');
    my $bp_end = $self->param_required('bp_end');
    my $overlap = $self->param('overlap');
    
    my $vcf = $self->param_required('input_vcf');
    my $depth_matrix_dir = $self->param_required('depth_matrix_dir');
    my $supp_dp_mx_dir = $self->param_required('supp_dp_mx_dir');
    my $sample_panel = $self->param_required('sample_panel');
    my $related_sample_list = $self->param_required('related_sample_list');
    my $simple_var_gl_file_dir = $self->param_required('simple_var_gl_file_dir');
    my $complex_var_gl_file_dir = $self->param_required('complex_var_gl_file_dir');
    
    my $tabix =  $self->param_required('tabix');   
	my $save_files_from_deletion = $self->param('save_files_from_deletion');

    my $slices = fai_to_slices(
          fai => $fai,
          SQ_start => $SQ_start, SQ_end => $SQ_end,
          bp_start => $bp_start, bp_end => $bp_end,
          );
    throw("expected exactly one slice, not ". scalar @$slices) if scalar @$slices != 1;

    $slices->[0]->extend($overlap);
   
    my $chrom = $slices->[0]->SQ_name,
    my $s;
    my $e;          
    if (!$slices->[0]->is_whole_SQ) {
      $s = $slices->[0]->start;
      $e = $slices->[0]->end;
    }
	
	my $chunk = $chrom . ":" . $s . "-" . $e;

	my $depth_matrix = $depth_matrix_dir . "/" . "ALL.chr" . $chrom . ".20140717.matrix.gz";
	my $supp_dp_mx = $supp_dp_mx_dir . "/" . "chr" . $chrom . ".20140718.matrix.gz";
	my $simple_var_gl_file = $simple_var_gl_file_dir . "/" . "ALL.chr" . $chrom . ".phase3_bc_union.20130502.biallelic_svm_snps_indelsAF0.005_svs.gl.vcf.gz";
	my $complex_var_gl_file = $complex_var_gl_file_dir . "/" . "ALL.chr" . $chrom . ".mvncall_sorted_off_by_one.20130502.non_biallelic_snps_svs_microsat.genotypes.vcf.gz";

	my $object = ReseqTrack::Tools::AddReadDpAddGlCalculateAF->new (
		-program					=> $tabix,
		-input_files				=> $vcf,
		-depth_matrix				=> $depth_matrix,
		-supp_dp_mx					=> $supp_dp_mx,
		-region						=> $chunk,
		-working_dir				=> $self->output_dir,
		-sample_panel				=> $sample_panel,
		-related_sample_list		=> $related_sample_list,
		-simple_var_gl_file			=> $simple_var_gl_file,
		-complex_var_gl_file		=> $complex_var_gl_file,
		-save_files_from_deletion	=> $save_files_from_deletion,
	);
	
	$self->run_program($object);

    $self->output_param('vcf' => $object->output_files->[0]);

}

1;