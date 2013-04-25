
package ReseqTrack::HiveProcess::RunVariantCall::CallByFreebayes;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::RunVariantCall::CallByFreebayes;
use ReseqTrack::HiveUtils::SequenceRegionsUtils qw(get_sequence_name);
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $bam = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $seq_index = $self->param('seq_index');
    my $region_start = $self->param('region_start');
    my $region_end = $self->param('region_end');
    throw("seq_index is an obligatory parameter") if ! defined $seq_index;
    throw("region_start is an obligatory parameter") if ! defined $region_start;
    throw("region_end is an obligatory parameter") if ! defined $region_end;

    my $sequence_name = get_sequence_name($self->data_dbc, $seq_index);

    $self->data_dbc->disconnect_when_inactive(1);

    my $variant_caller = ReseqTrack::Tools::RunVariantCall::CallByFreebayes->new(
          -input_files => $bam,
          -program => $self->param('freebayes'),
          -bgzip => $self->param('bgzip'),
          -working_dir => $self->output_dir,
          -job_name => $self->job_name,
          -chrom => $sequence_name,
          -region => "$region_start-$region_end",
          -reference => $self->param('reference'),
          -options => $self->param('options'),
          );
    $variant_caller->run;

    $self->output_this_branch('vcf' => $variant_caller->output_files);
    $self->data_dbc->disconnect_when_inactive(0);

}

1;

