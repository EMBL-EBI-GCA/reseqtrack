
package ReseqTrack::HiveProcess::SplitFastq;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);
use ReseqTrack::Tools::FileUtils qw(get_count_stats);
use ReseqTrack::Tools::RunSplit;
use POSIX qw(ceil);
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    $self->param_required('max_bases');
    $self->param_required('fastq');

    my $program_file = $self->param_is_defined('program_file') ? $self->param('program_file') || undef;

    my $output_dir = $self->output_dir;


    my $fastq_ids = $self->get_param_values('fastq');
    my $fastqs = $self->get_files($fastq_ids);
    my ($mate1, $mate2, $frag) = assign_files($fastqs);
    my ($mate1_path, $mate2_path, $frag_path) = map {$_ ? $_->name : ''} ($mate1, $mate2, $frag);
    throw ("No mate for $mate1_path") if ($mate1_path && ! $mate2_path);
    throw ("No mate for $mate2_path") if ($mate2_path && ! $mate1_path);

    my %line_count_hash;
    if ($mate1 && $mate2) {
      my ($read_count1, $base_count1) = get_count_stats($mate1);
      my ($read_count2, $base_count2) = get_count_stats($mate2);
      throw("read counts don't match: $mate1_path $read_count1 $mate2_path $read_count2")
            if ($read_count1 != $read_count2);

      my $num_output_files = ceil( 0.5 * ($base_count1 +$base_count2) / $max_bases );
      my $line_count = 4 * ceil($read_count1 / $num_output_files);
      $line_count_hash{$mate1_path} = $line_count;
      $line_count_hash{$mate2_path} = $line_count;
    }
    if ($frag) {
      my ($read_count, $base_count) = get_count_stats($frag);
      my $num_output_files = ceil($base_count / $max_bases);
      my $line_count = 4 * ceil($read_count / $num_output_files);
      $line_count_hash{$frag_path} = $line_count;
    }

    $self->data_dbc->disconnect_when_inactive(1);

    my $run_split = ReseqTrack::Tools::RunSplit->new(
        -program => $program_file,
        -working_dir => $output_dir,
        -line_count_hash => \%line_count_hash,
        );
    $run_split->run;
    my $output_file_hash = $run_split->output_file_hash;

    if ($mate1 && $mate2) {
      my @labels_mate1 = keys %{$output_file_hash->{$mate1_path}};
      foreach my $label (sort {$a <=> $b} @labels_mate1) {
        throw("no matching mate2 file with label $label") if (!$output_file_hash->{$mate2_path}->{$label});
        my @files = map {$output_file_hash->{$_}->{$label}} ($mate1_path, $mate2_path);
        my $file_ids = $self->make_file_ids(\@files);
        $self->prepare_child_output_id("$run_id.mate_chunk$label", {'fastq' => $file_ids});
      }
    }
    if ($frag) {
      foreach my $label (sort {$a <=> $b} keys %{$output_file_hash->{$frag_path}}) {
        my $file_path = $output_file_hash->{$frag_path}{$label};
        my $file_ids = $self->make_file_ids([$file_path]);
        $self->prepare_child_output_id("$run_id.frag_chunk$label", {'fastq' => $file_ids->[0]});
      }
    }

    $self->data_dbc->disconnect_when_inactive(0);

}

1;

