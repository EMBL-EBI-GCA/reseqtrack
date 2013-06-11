
package ReseqTrack::Hive::Process::RunVariantCall;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    $self->param_required('bam');
    my $bam = $self->get_param_values('bam');

    my $module_name = $self->param_required('module_name');
    my $fai = $self->param_required('fai');
    my $SQ_start = $self->param_required('SQ_start');
    my $SQ_end = $self->param_required('SQ_end');
    my $bp_start = $self->param_required('bp_start');
    my $bp_end = $self->param_required('bp_end');
    my $overlap = $self->param_is_defined('overlap') ? $self->param('overlap') : 0;

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
      $module_args{'-program'} = $self->param_is_defined('samtools') ? $self->param('samtools') : undef;
      $module_args{'-bcftools'} = $self->param_is_defined('bcftools') ? $self->param('bcftools') : undef;
      $module_args{'-vcfutils'} = $self->param_is_defined('vcfutils') ? $self->param('vcfutils') : undef;
      $module_args{'-bgzip'} = $self->param_is_defined('bgzip') ? $self->param('bgzip') : undef;
    }
    elsif ($module_name eq 'CallByGATK') {
      $module_args{'-java_exe'} = $self->param_is_defined('java_exe') ? $self->param('java_exe') : undef;
      $module_args{'-jvm_args'} = $self->param_is_defined('jvm_args') ? $self->param('jvm_args') : undef;
      $module_args{'-gatk_path'} = $self->param_is_defined('gatk_dir') ? $self->param('gatk_dir') : undef;
    }
    elsif ($module_name eq 'CallByFreebayes') {
      $module_args{'-program'} = $self->param_is_defined('freebayes') ? $self->param('freebayes') : undef;
    }


    my $variant_caller = $module->new(
          -input_files => $bam,
          -working_dir => $self->output_dir,
          -job_name => $self->job_name,
          -chrom => $slices->[0]->SQ_name,
          -reference => $self->param_required('reference'),
          -options => $self->param_is_defined('options') ? $self->param('options') : undef,
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

