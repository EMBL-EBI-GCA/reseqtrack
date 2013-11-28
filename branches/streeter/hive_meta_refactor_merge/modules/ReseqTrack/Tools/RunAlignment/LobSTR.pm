package ReseqTrack::Tools::RunAlignment::LobSTR;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);


use base qw(ReseqTrack::Tools::RunAlignment);

sub DEFAULT_OPTIONS { return {
        'threads' => 1,
        };
}

sub new {

    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my ( $ref_index_prefix)
        = rearrange( [
            qw( REF_INDEX_PREFIX
                )], @args);


    $self->ref_index_prefix($ref_index_prefix);


    return $self;

}
#################################################################


sub run_alignment {
    my ($self) = @_;

    throw("need a ref_index_prefix") if ! $self->ref_index_prefix;

    TYPE:
    foreach my $aln_type ('MATE', 'FRAG') {
      next TYPE if $aln_type eq 'MATE' && (!$self->mate1_file || !$self->mate2_file);
      next TYPE if $aln_type eq 'FRAG' && !$self->fragment_file;

      my $output_prefix = $self->working_dir() . '/'
          . $self->job_name . '_se';
      $output_prefix = $self->working_dir() . '/' . $self->job_name;
      $output_prefix .= ($aln_type eq 'MATE' ? '_pe' : '_se');

      my $output_file = "$output_prefix.aligned.bam";
      $output_file =~ s{//}{/};

      my @cmd_words;
      push(@cmd_words, $self->program);

      if ($aln_type eq 'MATE') {
        my $fastq_string_1 = $self->get_fastq_cmd_string('mate1');
        my $fastq_string_2 = $self->get_fastq_cmd_string('mate2');
        push(@cmd_words, '--p1', $fastq_string_1, '--p2', $fastq_string_2);
        if ($fastq_string_1 =~ /\.gz\S*$/) {
          push(@cmd_words, '--gzip');
        }
      }
      else {
        my $fastq_string = $self->get_fastq_cmd_string('frag');
        push(@cmd_words, '--files', $fastq_string);
        if ($fastq_string =~ /\.gz\S*$/) {
          push(@cmd_words, '--gzip');
        }
      }

      push(@cmd_words, '--index-prefix', $self->ref_index_prefix);
      push(@cmd_words, '--threads', $self->options('threads') || 1);
      push(@cmd_words, '--fastq');

      if (my $sample = $self->read_group_fields->{'SM'}) {
        push(@cmd_words, '--rg-sample', $sample);
      }
      if (my $library = $self->read_group_fields->{'LB'}) {
        push(@cmd_words, '--rg-lib', $library);
      }

      push(@cmd_words, '--out', $output_prefix);

      my $cmd_line = join(' ', @cmd_words);

      $self->output_files($output_file);
      $self->created_files("$output_prefix.stats");
      $self->execute_command_line($cmd_line);
    }

}


sub ref_index_prefix {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{ref_index_prefix} = $arg;
    }
    return $self->{ref_index_prefix};
}

1;

