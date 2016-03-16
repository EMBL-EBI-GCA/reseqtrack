package ReseqTrack::Hive::Process::RunBismark;

use strict;
use warnings;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::RunBismark;
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

    my $command = $self->param_required('command');

    my @allowed_cmds = qw(aln methext);
    throw("Don't recognise command $command. Acceptable commands are: @allowed_cmds")
	if (! grep {$command eq $_ } @allowed_cmds);


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

    if ($command eq 'aln') {
	my $chunk_label = $self->param_required('chunk_label');
	my $reference = $self->param_required('reference');
	my $fastqs = $self->param_as_array('fastq');
	foreach my $fastq (@$fastqs) {
	    check_file_exists($fastq);
	}
	my $run_alignment;
	if (scalar(@$fastqs)>1) {
	    $run_alignment = ReseqTrack::Tools::RunBismark->new(
		-base => $chunk_label,
                -mate1_file => $fastqs->[0],
		-mate2_file => $fastqs->[1],
                -program => $self->param('program_file'),
                -samtools => $self->param('samtools'),
                -output_format => 'BAM',
                -working_dir => $self->output_dir,
                -reference => $reference,
                -job_name => $self->job_name,
                -options => $self->param('options'),
                );
	} else {
	    $run_alignment = ReseqTrack::Tools::RunBismark->new(
		-base => $chunk_label,
		-fragment_file => $fastqs->[0],
		-program => $self->param('program_file'),
		-samtools => $self->param('samtools'),
		-output_format => 'BAM',
		-working_dir => $self->output_dir,
		-reference => $reference,
		-job_name => $self->job_name,
		-options => $self->param('options'),
		);
	}

	$self->run_program($run_alignment,$command);

	$self->output_param('bam', $run_alignment->output_files->[0]);
    } elsif ($command eq 'methext') {
	my $fastqs = $self->param_as_array('fastq');
	my $runmode;
	if (scalar($fastqs)>1) {
	    $runmode='PAIRED';
	} else {
	    $runmode='SINGLE';
	}
	my $bamfile = $self->param_required('bam');
	check_file_exists($bamfile);
	my $run_methext=ReseqTrack::Tools::RunBismark->new(
	    -input_files => $bamfile,
	    -working_dir => $self->output_dir,
	    -runmode => $runmode
	    );
	$self->run_program($run_methext,$command);

    }
}

1;

