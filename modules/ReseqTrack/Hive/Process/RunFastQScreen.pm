
package ReseqTrack::Hive::Process::RunFastQScreen;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
#use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
#use ReseqTrack::Tools::AttributeUtils qw( create_attribute_for_object );
use ReseqTrack::Tools::QC::FastQScreen;


sub param_defaults {
  return {
    store_attributes => 0,
    program_file => undef,
  };
}


sub run {
    my $self = shift @_;

    $self->param_required('fastq');

    my $fastqs = $self->param_as_array('fastq');
    throw("Expecting one fastq file") if scalar @$fastqs != 1;

    my $output_dir = $self->output_dir;
    my $conf_file = $self->param('conf_file');
    my $fastqscreen = ReseqTrack::Tools::QC::FastQScreen->new(
      -conf_file => $conf_file,
      -job_name => $self->job_name,
      -program => $self->param('program_file'),
      -keep_text => 1,
      -keep_summary => 1,
      -keep_zip => 1,
      -working_dir => $output_dir,
      -input_files => $fastqs,
    );

    $self->run_program($fastqscreen);
}

1;
