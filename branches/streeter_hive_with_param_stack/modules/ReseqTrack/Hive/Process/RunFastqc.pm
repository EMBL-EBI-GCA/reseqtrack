
package ReseqTrack::Hive::Process::RunFastqc;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::AttributeUtils qw( create_attribute_for_object );
use ReseqTrack::Tools::QC::FastQC;


sub param_defaults {
  return {
    store_attributes => 0,
    program_file => undef,
  };
}


sub run {
    my $self = shift @_;

    $self->param_required('fastq');
    my $store_attributes = $self->param('store_attributes') ? 1 : 0;
    my $db_params = $store_attributes ? $self->param_required('reseqtrack_db') : undef;

    my $fastqs = $self->file_param_to_flat_array('fastq');
    throw("Expecting one fastq file") if scalar @$fastqs != 1;
    my $fastq = $fastqs->[0];

    my ($db, $fastq_object);
    if ($store_attributes) {
      $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});
      $fastq_object = $db->get_FileAdaptor->fetch_by_name($fastq);
      throw("did not get file with name $fastq") if !$fastq_object;
      $db->dbc->disconnect_when_inactive(1);
    }

    my $output_dir = $self->output_dir;
    my $fastqc = ReseqTrack::Tools::QC::FastQC->new(
      -job_name => $self->job_name,
      -program => $self->param('program_file'),
      -keep_text => 1,
      -keep_summary => 1,
      -keep_zip => 1,
      -working_dir => $output_dir,
      -input_files => $fastqs,
    );

    $self->run_program($fastqc);

    my @summary_paths = grep { /_summary\.txt$/ } @{$fastqc->output_files};
    my @report_paths = grep { /_report\.txt$/ } @{$fastqc->output_files};
    my @zip_paths = grep { /\.zip$/ } @{$fastqc->output_files};

    throw("unexpected number of summary paths: " . scalar @summary_paths) if @summary_paths != 1;
    throw("unexpected number of report paths: " . scalar @report_paths) if @report_paths != 1;
    throw("unexpected number of zip paths: " . scalar @zip_paths) if @zip_paths != 1;
    my ($summary_path, $report_path, $zip_path) = ($summary_paths[0], $report_paths[0], $zip_paths[0]);

    if ($store_attributes) {
      $db->dbc->disconnect_when_inactive(0);
      my $attributes = $fastq_object->attributes;
      
      open (my $summary_fh, '<', $summary_path) or throw("Could not open $summary_path $!");
      while (<$summary_fh>){
          chomp;
          my ($value,$key,$name) = split /\t/;
          push @$attributes, create_attribute_for_object($fastq_object,"FASTQC:$key",$value) if ($value && $key);
      }
      
      $fastq_object->uniquify_attributes($attributes);
      $db->get_FileAdaptor->store_attributes($fastq_object);
    }

    $self->output_param('fastqc_summary', $summary_path);
    $self->output_param('fastqc_report', $report_path);
    $self->output_param('fastqc_zip', $zip_path);
}

1;

