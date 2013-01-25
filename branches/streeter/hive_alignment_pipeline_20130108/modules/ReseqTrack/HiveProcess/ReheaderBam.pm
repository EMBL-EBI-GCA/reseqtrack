
package ReseqTrack::HiveProcess::ReheaderBam;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::ReheaderBam;
use File::Basename qw(fileparse);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my ($self) = @_;

    my $bams = $self->param('bam') || die "'bam' is an obligatory parameter";
    my $fastqs = $self->param('fastq') || die "'fastqs' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";
    my $job_name = $self->param('job_name') or die "'job_name' is an obligatory parameter";
    my $samtools = $self->param('samtools');

    my @extra_header_lines;
    $fastqs = ref($fastqs) eq 'ARRAY' ? $fastqs : [$fastqs];
    foreach my $fastq (grep {$_} @$fastqs) {
      my $basename = fileparse($fastq);
      push(@extra_header_lines, "\@CO\tFASTQ=$basename");
    }

    my $reheader_object = ReseqTrack::Tools::ReheaderBam->new(
      -input_files  => $bams,
      -working_dir  => $output_dir,
      -job_name     => $job_name,
      -samtools     => $samtools,
      -header_lines_file => $self->param('header_lines_file'),
      -extra_header_lines => \@extra_header_lines,
      -options      => {reuse_old_header => 1,
                        replace_PG => 1,
                        replace_CO => 1},
      -SQ_fields_hash => {AS => $self->param('SQ_assembly'),
                          SP => $self->param('SQ_species'),
                          UR => $self->param('SQ_uri')},
    );
    $reheader_object->run;
    $self->output_this_branch('bam' => $reheader_object->output_files);
}


1;

