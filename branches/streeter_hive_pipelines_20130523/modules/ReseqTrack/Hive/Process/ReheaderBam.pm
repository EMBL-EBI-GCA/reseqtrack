
package ReseqTrack::Hive::Process::ReheaderBam;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::ReheaderBam;
use File::Basename qw(fileparse);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my ($self) = @_;

    $self->param_required('bam');
    my $bams = $self->file_param_to_flat_array('bam');
    $self->param_required('fastq');
    my $fastqs = $self->file_param_to_flat_array('fastq');

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
      -program      => $self->param_is_defined('samtools') ? $self->param('samtools') : undef,
      -header_lines_file => $self->param_is_defined('header_lines_file') ? $self->param('header_lines_file') : undef,
      -dict_file => $dict_file,
      -extra_header_lines => \@extra_header_lines,
      -options      => {reuse_old_header => 1,
                        replace_PG => 1,
                        replace_CO => 1},
      -SQ_fields_hash => {AS => $self->param_is_defined('SQ_assembly') ? $self->param('SQ_assembly') : undef,
                          SP => $self->param_is_defined('SQ_species') ? $self->param('SQ_species') : undef,
                          UR => $self->param_is_defined('SQ_uri') ? $self->param('SQ_uri') : undef},
    );

    $self->run_program($reheader_object);

    my $output_files = ref($self->param('bam')) ? $reheader_object->output_files : $reheader_object->output_files->[0];
    $self->output_param('bam' => $output_files);
}


1;

