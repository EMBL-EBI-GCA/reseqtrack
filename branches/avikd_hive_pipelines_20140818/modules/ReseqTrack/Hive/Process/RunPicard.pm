
package ReseqTrack::Hive::Process::RunPicard;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::RunPicard;
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use Data::Dump qw(dump); ## fix

sub param_defaults {
  return {
    java_exe => undef,
    jvm_args => undef,
    picard_dir => undef,
    create_index => undef,
    options => {},
    add_attributes  => 0,
  };
}


sub run {
    my $self = shift @_;
    $self->param_required('bam');
    my $bams = $self->param_as_array('bam');
    my $command = $self->param_required('command');
    
    my $add_attributes = $self->param( 'add_attributes' ) ? 1 : 0;  

    my @allowed_cmds = ReseqTrack::Tools::RunPicard->get_valid_commands;
    throw( "Don't recognise command $command. Acceptable commands are: @allowed_cmds")
      if ( !grep { $command eq $_ } @allowed_cmds );

    my %options = (validation_stringency => 'SILENT');
    my $bam_name_length = 0;
    foreach my $bam_name (@$bams) {
      $bam_name_length += length($bam_name);
    }
    if ($bam_name_length > 120000) {
      $options{'shorten_input_names'} => 1;
    }
    %options = (%options, %{$self->param('options')});

    my $picard_object = ReseqTrack::Tools::RunPicard->new(
      -input_files  => $bams,
      -working_dir  => $self->output_dir,
      -job_name     => $self->job_name,
      -java_exe     => $self->param('java_exe'),
      -jvm_options  => $self->param('jvm_args'),
      -picard_dir     => $self->param('picard_dir'),
      -options      => \%options,
      -create_index => $self->param('create_index'),
      -keep_metrics => 0,
    );

    $self->run_program($picard_object, $command);

    my $output_bams = $picard_object->output_bam_files;
    my $output_bais = $picard_object->output_bai_files;

    $self->output_param('bam', $output_bams);
    $self->output_param('bai', $output_bais);
   
    if( $add_attributes && $command eq 'mark_duplicates'){  ## now only for mark_duplicates command
      my $generated_metrics = $picard_object->output_metrics_object;
      throw('metrics object not found') unless $generated_metrics;
      dump( $generated_metrics );
      $self->output_param('attribute_metrics', $generated_metrics); 
    }
}


1;

