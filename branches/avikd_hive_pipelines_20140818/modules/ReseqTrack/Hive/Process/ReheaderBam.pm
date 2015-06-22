
package ReseqTrack::Hive::Process::ReheaderBam;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::ReheaderBam;
use File::Basename qw(fileparse);


sub param_defaults {
  return {
    samtools => undef,
    header_lines_file => undef,
    SQ_assembly => undef,
    SQ_species => undef,
    SQ_uri => undef,
  };
}

sub run {
    my ($self) = @_;

    $self->param_required('bam');
    my $bams = $self->param_as_array('bam');
    $self->param_required('fastq');
    my $fastqs = $self->param_as_array('fastq');

    my @extra_header_lines;
    foreach my $fastq (grep {$_} @$fastqs) {
      my $basename = fileparse($fastq);
      push(@extra_header_lines, "\@CO\tFASTQ=$basename");
    }

    my $dict_file;
    if ($self->param_is_defined('dict_file')) {
      $dict_file = $self->param('dict_file');
    }
    elsif ($self->param_is_defined('reference')) {
      $dict_file = $self->param('reference');
      $dict_file =~ s/fa(?:sta)?(?:\.gz)?$/dict/;
    }

    my $reheader_object = ReseqTrack::Tools::ReheaderBam->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -job_name     => $self->job_name,
      -program      => $self->param('samtools'),
      -header_lines_file => $self->param('header_lines_file'),
      -dict_file => $dict_file,
      -extra_header_lines => \@extra_header_lines,
      -options      => {reuse_old_header => 1,
                        replace_PG => 1,
                        replace_CO => 1},
      -SQ_fields_hash => {AS => $self->param('SQ_assembly'),
                          SP => $self->param('SQ_species'),
                          UR => $self->param('SQ_uri')},
    );

    $self->run_program($reheader_object);

    my $output_files = ref($self->param('bam')) ? $reheader_object->output_files : $reheader_object->output_files->[0];
    $self->output_param('bam' => $output_files);
}


1;

