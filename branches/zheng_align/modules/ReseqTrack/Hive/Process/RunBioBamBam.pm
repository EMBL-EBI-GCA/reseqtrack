package ReseqTrack::Hive::Process::RunBioBamBam;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunBioBamBam;
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);

sub param_defaults {
  return {
    biobambam_dir		=> undef,
    create_index		=> undef,
    remove_duplicates	=> undef,
    sort_order			=> undef,
    options				=> {},
  };
}

sub run {
    my $self = shift @_;
    
    $self->param_required('bam');
    my $bams = $self->param_as_array('bam');
    my $command = $self->param_required('command');

    my %options;
    
    %options = (%options, %{$self->param('options')});

    my $biobambam_object = ReseqTrack::Tools::RunBioBamBam->new(
      -input_files			=> $bams,
      -working_dir			=> $self->output_dir,
      -job_name				=> $self->job_name,
      -biobambam_dir		=> $self->param('biobambam_dir'),
      -options				=> \%options,
      -create_index			=> $self->param('create_index'),
      -sort_order			=> $self->param('sort_order'),
      -command				=> $command,
      -keep_metrics			=> 0,
      -remove_duplicates	=> 0.
    );

    $self->run_program($biobambam_object);

    my $output_bams = $biobambam_object->output_bam_files;
    my $output_bais = $biobambam_object->output_bai_files;

    $self->output_param('bam', $output_bams);
    $self->output_param('bai', $output_bais);

}


1;

