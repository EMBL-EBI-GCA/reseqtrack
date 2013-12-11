
package ReseqTrack::Hive::Process::BWA;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::RunAlignment::BWA;
use ReseqTrack::Tools::Exception qw(throw);

sub param_defaults {
  return {
    program_file => undef,
    samtools => undef,
    options => {},

    RGID => undef,
    RGCN => undef,
    RGLB => undef,
    RGPI => undef,
    RGSM => undef,
    RGDS => undef,
    RGPU => undef,
    RGPL => undef,
    sample_source_id => undef,
    center_name => undef,
    library_name => undef,
    paired_nominal_length => undef,
    run_source_id => undef,
    instrument_platform => undef,
  };
}

sub run {
    my $self = shift @_;

    $self->param_required('fastq');
    my $run_id = $self->param_required('run_id');
    my $reference = $self->param_required('reference');

    my $fastqs = $self->param_as_array('fastq');

    my %read_group_fields = (
      ID => $self->param('RGID') // $self->param('run_source_id'),
      CN => $self->param('RGCN') // $self->param('center_name'),
      LB => $self->param('RGLB') // $self->param('library_name'),
      PI => $self->param('RGPI') // $self->param('paired_nominal_length'),
      SM => $self->param('RGSM') // $self->param('sample_alias'),
      DS => $self->param('RGDS') // $self->param('study_source_id'),
      PU => $self->param('RGPU') // $self->param('run_source_id'),
      PL => $self->param('RGPL') // $self->param('instrument_platform'),
    );

    foreach my $fastq (@$fastqs) {
      check_file_exists($fastq);
    }


    my $run_alignment = ReseqTrack::Tools::RunAlignment::BWA->new(
          -input_files => $fastqs,
          -program => $self->param('program_file'),
          -samtools => $self->param('samtools'),
          -output_format => 'BAM',
          -working_dir => $self->output_dir,
          -reference => $reference,
          -job_name => $self->job_name,
          -paired_length => $read_group_fields{'PI'},
          -read_group_fields => \%read_group_fields,
          -options => $self->param('options'),
          );

    $self->run_program($run_alignment);

    $self->output_param('bam', $run_alignment->output_files->[0]);

}

1;

