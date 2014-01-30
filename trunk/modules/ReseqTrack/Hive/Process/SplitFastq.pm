
package ReseqTrack::Hive::Process::SplitFastq;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);
use ReseqTrack::Tools::FileUtils qw(get_count_stats);
use ReseqTrack::Tools::RunSplit;
use POSIX qw(ceil);
use ReseqTrack::Tools::Exception qw(throw);

sub param_defaults {
  return {
    program_file => undef,
    samtools => undef,
    run_alias => undef,
  };
}


sub run {
    my $self = shift @_;

    my $max_reads = $self->param_required('max_reads');
    $self->param_required('fastq');

    my $program_file = $self->param('program_file');

    my $output_dir = $self->output_dir;

    my $run_source_id = $self->param_required('run_source_id');
    my $run_alias = $self->param('run_alias');
    my $search_string = ($run_source_id && $run_alias) ? "$run_source_id|$run_alias" : $run_source_id || $run_alias;
    my @regexs = (qr/(?:$search_string})_1\.(?:\w+\.)*f(?:ast)?q(?:\.gz)?/i,
                  qr/(?:$search_string})_2\.(?:\w+\.)*f(?:astq)?(?:\.gz)?/i,
                  qr/(?:$search_string})\.(?:\w+\.)*f(?:ast)?q(?:\.gz)?/i);


    my $fastqs = $self->param_as_array('fastq');
    my ($mate1, $mate2, $frag) = assign_files($fastqs, \@regexs);
    throw ("No mate for $mate1") if ($mate1 && ! $mate2);
    throw ("No mate for $mate2") if ($mate2 && ! $mate1);

    my $run_split = ReseqTrack::Tools::RunSplit->new(
        -input_files => $fastqs,
        -program => $program_file,
        -working_dir => $output_dir,
        -line_count => 4*$max_reads
        );
    $self->run_program($run_split);
    my $output_file_hash = $run_split->output_file_hash;

    if ($mate1 && $mate2) {
      my @labels_mate1 = keys %{$output_file_hash->{$mate1}};
      foreach my $label (sort {$a <=> $b} @labels_mate1) {
        throw("no matching mate2 file with label $label") if (!$output_file_hash->{$mate2}->{$label});
        my @files = map {$output_file_hash->{$_}->{$label}} ($mate1, $mate2);
        $self->prepare_factory_output_id({'fastq' => \@files, 'chunk' => "mate_chunk$label"});
      }
    }
    if ($frag) {
      foreach my $label (sort {$a <=> $b} keys %{$output_file_hash->{$frag}}) {
        my $file_path = $output_file_hash->{$frag}{$label};
        $self->prepare_factory_output_id({'fastq' => $file_path, 'chunk' => "frag_chunk$label"});
      }
    }

}

1;

