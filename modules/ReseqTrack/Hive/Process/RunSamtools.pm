
package ReseqTrack::Hive::Process::RunSamtools;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use ReseqTrack::Tools::RunSamtools;
use File::Basename;

sub param_defaults {
  return {
    samtools_options => {},
    program_file => undef,
    reference => undef,
  };
}


sub run {
    my $self = shift @_;
    $self->param_required('bam');
    my $bams = $self->param_as_array('bam');
    my $command = $self->param_required('command');

    my @allowed_cmds = qw(merge sort index fix_and_calmd calmd fixmate sam_to_bam bam_to_cram);
    throw("Don't recognise command $command. Acceptable commands are: @allowed_cmds")
      if (! grep {$command eq $_ } @allowed_cmds);

    my $options = $self->param('samtools_options');
    my @allowed_options = keys %{&ReseqTrack::Tools::RunSamtools::DEFAULT_OPTIONS};
    foreach my $option (keys %$options) {
      throw("Don't recognise option $option. Acceptable options are: @allowed_options")
        if (! grep {$option eq $_ } @allowed_options);
    }


    my $samtools_object = ReseqTrack::Tools::RunSamtools->new(
      -input_files  => $bams,
      -program      => $self->param('program_file'),
      -working_dir  => $self->output_dir,
      -job_name     => $self->job_name,
      -reference    => $self->param('reference'),
      -options      => $options,
    );

    $self->run_program($samtools_object, $command);

    my $output_files = $samtools_object->output_files;
	
	if (basename($output_files->[0]) =~ /\.cram/i  || basename($output_files->[0]) =~ /\.crai/i) {
		print "cram output is $output_files->[0]\n";
		$self->output_param($command eq 'index' ? 'crai' : 'cram'  => $output_files);
	}
	else {
		print "bam output is $output_files->[0]\n";
	    $self->output_param($command eq 'index' ? 'bai' : 'bam'  => $output_files);
	}
		
}


1;

