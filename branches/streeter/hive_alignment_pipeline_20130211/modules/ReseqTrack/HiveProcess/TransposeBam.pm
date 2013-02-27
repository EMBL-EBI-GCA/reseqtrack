
package ReseqTrack::HiveProcess::TransposeBam;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::RunTransposeBam;
use ReseqTrack::HiveUtils::SequenceRegionsUtils qw(get_region_strings);
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $bams = $self->param('bams') || throw("'bam' is an obligatory parameter");
    my $seq_index_start = $self->param('seq_index_start');
    my $seq_index_end = $self->param('seq_index_end');
    throw("'seq_index_start' is an obligatory parameter") if !defined $seq_index_start;
    throw("'seq_index_end' is an obligatory parameter") if !defined $seq_index_end;

    my $regions = get_region_strings($self->data_dbc, $seq_index_start, $seq_index_end);


    $self->data_dbc->disconnect_when_inactive(1);

    my $bam_transposer = ReseqTrack::Tools::RunTransposeBam->new(
          -input_files => $bams,
          -program => $self->param('program_file'),
          -working_dir => $self->output_dir,
          -job_name => $self->job_name,
          -regions => $regions,
          -options => {'build_index' => $self->param('create_index')},
          );
    $bam_transposer->run;

    my %output_data;
    $output_data{'bam'} = [grep { /\.bam$/ } @{$bam_transposer->output_files}];
    if ($self->param('create_index')) {
      $output_data{'bai'} = [grep { /\.bai$/ } @{$bam_transposer->output_files}];
    }

    $self->output_this_branch(%output_data);
    $self->data_dbc->disconnect_when_inactive(0);

}

1;

