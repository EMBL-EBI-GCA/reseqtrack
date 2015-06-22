
package ReseqTrack::Hive::Process::RunRnaAlignment;

use strict;

use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);


sub param_defaults {
  return {
    options => undef,
    samtools => undef,
    gsnap => undef,
    star => undef,
  };
}


sub run {
    my $self = shift @_;
    
    $self->param_required('fastq');
    my $run_id = $self->param_required('run_id');
    my $reference = $self->param_required('reference');

    my $fastqs = $self->param_as_array('fastq');
    
    my $module_name = $self->param_required('module_name');
   
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

    foreach my $fastq (@$fastqs) {
      check_file_exists($fastq);
    }

    my $module = load_module($module_name);

    my %module_args;
    if ($module_name eq 'GSNAP') {
      $module_args{'-program'} = $self->param('gsnap');
    }
    elsif ($module_name eq 'STAR') {
      $module_args{'-program'} = $self->param('star');
    }

    my $run_alignment = $module->new(
          -input_files => $fastqs,
          -program => $self->param('program_file'),
          -working_dir => $self->output_dir,
          -output_format => 'BAM',
          -job_name => $self->job_name,
          -reference => $self->param_required('reference'),
          -read_group_fields => \%read_group_fields,
          -samtools => $self->param('samtools'),
          -options => $self->param('options'),
          %module_args,
          );
    

    $self->run_program($run_alignment);
    $self->output_param('bam', $run_alignment->output_bam_files->[0] );
    
}

sub load_module {
  my ($module_name) = @_;
  my $file = "ReseqTrack/Tools/RunAlignment/$module_name.pm";
  eval {
    require "$file";
  };
  if ($@) {
    throw("cannot load $file: $@")
  }
  my $module = "ReseqTrack::Tools::RunAlignment::$module_name";
  return $module;
}

1;

