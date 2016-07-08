
package AccessibleGenome::BamUtilHiveProcess;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);

use AccessibleGenome::RunBamUtil;


sub param_defaults {
  return {
    options => undef,

    bam_util_exe => '/nfs/1000g-work/G1K/work/bin/bamUtil/bin/bam',
  };
}


sub run {
    my $self = shift @_;

    my $bam = $self->param_required('bam');
    my $fai = $self->param_required('fai');

    my $SQ= $self->param_required('SQ');
    my $bp_start = $self->param_required('bp_start');
    my $bp_end = $self->param_required('bp_end');

    my $bam_util = AccessibleGenome::RunBamUtil->new(
          -input_files => $bam,
          -working_dir => $self->output_dir,
          -program => $self->param('bam_util_exe') // $self->param_defaults->{'bam_util_exe'},
          -job_name => $self->job_name,
          -SQ => $SQ,
          -region_start => $bp_start,
          -region_end => $bp_end,
          -options => $self->param('options'),
          );

    $self->run_program($bam_util);

    $self->output_param('stats' => $bam_util->output_files->[0]);

}

1;

