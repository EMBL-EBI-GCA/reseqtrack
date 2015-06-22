
package ReseqTrack::Hive::Process::RunVep;

use strict;
use ReseqTrack::Tools::RunVep;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw(fai_to_slices bed_to_slices);
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    $self->param_required('vcf');
    my $vcf = $self->file_param_to_flat_array('vcf');
    my $fai = $self->param_required('fai');
    
    $self->param_required('bed');
    my $beds = $self->file_param_to_flat_array('bed');
    my $bed = $$beds[0];
    
    my $SQ_start = $self->param_required('SQ_start');
    my $SQ_end = $self->param_required('SQ_end');
    my $bp_start = $self->param_required('bp_start');
    my $bp_end = $self->param_required('bp_end');
    my $overlap = $self->param_is_defined('region_overlap') ? $self->param('region_overlap') : 0;
    my $create_index = $self->param_is_defined('create_index') ? $self->param('create_index') : 1;

    my $slices = fai_to_slices(
          fai => $fai,
          SQ_start => $SQ_start, SQ_end => $SQ_end,
          bp_start => $bp_start, bp_end => $bp_end,
          );
    if (defined $bed) {
      $slices = bed_to_slices(bed => $bed, parent_slices => $slices);
    }

    foreach my $slice (@$slices) {
      $slice->extend($overlap);
    }
    #$slices = join_overlapping_slices(slices => $slices, separation => 500);

    my @regions;
    foreach my $slice (@$slices) {
      my $region = $slice->SQ_name;
      if (!$slice->is_whole_SQ) {
        $region .= ':' . $slice->start . '-' . $slice->end;
      }
      push(@regions, $region);
    }
#########
throw("expected exactly one slice, not ". scalar @$slices) if scalar @$slices != 1;

my $vep_runner = ReseqTrack::Tools::RunVep->new(
      -input_files  => $vcf,
      -working_dir  => $self->output_dir,
      -program      => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
      -region => $slices->[0]->SQ_name . ":" . $slices->[0]->start . "-" . $slices->[0]->end,
      -options => $self->param_is_defined('options') ? $self->param('options') : undef,
      -tabix => $self->param_is_defined('tabix_exe') ? $self->param('tabix_exe') : undef,
      
    );

    $self->run_program($vep_runner);

    $self->output_param('vcf' => $vep_runner->output_files->[0]);

}

1;

