
package ReseqTrack::Hive::Process::RunVariantCall;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);
use ReseqTrack::Tools::Exception qw(throw);


sub param_defaults {
  return {
    overlap => 0,
    options => undef,

    samtools => undef,
    bcftools => undef,
    vcfutils => undef,
    bgzip => undef,

    java_exe => undef,
    jvm_args => undef,
    gatk_dir => undef,

    freebayes => undef,

    lobstr => undef,
  };
}


sub run {
    my $self = shift @_;

    $self->param_required('bam');
    my $bam = $self->param_as_array('bam');

    my $module_name = $self->param_required('module_name');
    my $fai = $self->param_required('fai');
    my $SQ_start = $self->param_required('SQ_start');
    my $SQ_end = $self->param_required('SQ_end');
    my $bp_start = $self->param_required('bp_start');
    my $bp_end = $self->param_required('bp_end');
    my $overlap = $self->param('overlap');

    my $module = load_module($module_name);

    my $slices = fai_to_slices(
          fai => $fai,
          SQ_start => $SQ_start, SQ_end => $SQ_end,
          bp_start => $bp_start, bp_end => $bp_end,
          );
    throw("expected exactly one slice, not ". scalar @$slices) if scalar @$slices != 1;

    $slices->[0]->extend($overlap);

    my %module_args;
    if ($module_name eq 'CallBySamtools') {
      $module_args{'-program'} = $self->param('samtools');
      $module_args{'-bcftools'} = $self->param('bcftools');
      $module_args{'-vcfutils'} = $self->param('vcfutils');
      $module_args{'-bgzip'} = $self->param('bgzip');
    }
    elsif ($module_name eq 'CallByBcftools') {
	$module_args{'-program'} = $self->param('bcftools');
	$module_args{'-samtools'} = $self->param('samtools');
    }
    elsif ($module_name eq 'CallByGATK') {
      $module_args{'-java_exe'} = $self->param('java_exe');
      $module_args{'-jvm_args'} = $self->param('jvm_args');
      $module_args{'-gatk_path'} = $self->param('gatk_dir');
    }
    elsif ($module_name eq 'CallByFreebayes') {
      $module_args{'-program'} = $self->param('freebayes');
      $module_args{'-bgzip'} = $self->param('bgzip');
    }
    elsif ($module_name eq 'CallByLobSTR') {
      $module_args{'-program'} = $self->param('lobstr');
      $module_args{'-noise_model'} = $self->param_required('noise_model');
      $module_args{'-str_info'} = $self->param_required('str_info');
      $module_args{'-ref_index_prefix'} = $self->param_required('ref_index_prefix');
      $module_args{'-bgzip'} = $self->param('bgzip');
      $self->param('reference', $self->param('ref_index_prefix'));
    }


    my $variant_caller = $module->new(
          -input_files => $bam,
          -working_dir => $self->output_dir,
          -job_name => $self->job_name,
          -chrom => $slices->[0]->SQ_name,
          -reference => $self->param_required('reference'),
          -options => $self->param('options'),
          %module_args,
          );
    if (!$slices->[0]->is_whole_SQ) {
      $variant_caller->region_start($slices->[0]->start);
      $variant_caller->region_end($slices->[0]->end);
    }

    $self->run_program($variant_caller);

    $self->output_param('vcf' => $variant_caller->output_files->[0]);

}

sub load_module {
  my ($module_name) = @_;
  my $file = "ReseqTrack/Tools/RunVariantCall/$module_name.pm";
  eval {
    require "$file";
  };
  if ($@) {
    throw("cannot load $file: $@")
  }
  my $module = "ReseqTrack::Tools::RunVariantCall::$module_name";
  return $module;
}

1;

