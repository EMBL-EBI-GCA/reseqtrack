=pod

=head1 NAME

ReseqTrack::Tools::RunPicard

=head1 SYNOPSIS

This is a class for running picard
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunPicard;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-picard_dir]   :
      string, directory containing picard jar files
  Arg [-program]   :
      string, the java executable
  Arg [-jvm_options]   :
      string, options for java, default is '-Xmx2g'
  Arg [-options]   :
      string, options for the picard jar file
  Arg [-tmp_dir]   :
      string, path of tmp_dir for picard, default is working_dir
  Arg [-index_bai]   :
      boolean, flag to produce an index file with suffix .bai
  Arg [-index_bam_bai]   :
      boolean, flag to produce an index file with suffix .bam.bai
  Arg [-replace_files]   :
      boolean, default 1, any input file will be deleted after output is created.
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunPicard object.
  Returntype: ReseqTrack::Tools::RunPicard
  Exceptions: 
  Example   : my $run_picard = ReseqTrack::Tools::RunPicard->new(
                -input_files => ['/path/sam1', '/path/sam2'],
                -picard_dir => '/path/to/picard/',
                -working_dir => '/path/to/dir/',
                -index_bam_bai => 1);

=cut


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $picard_dir, $jvm_options,
        $index_bai, $index_bam_bai,
        $options, $tmp_dir,
        $replace_files)
    = rearrange( [
         qw( PICARD_DIR JVM_OPTIONS
                INDEX_BAI INDEX_BAM_BAI
                OPTIONS TMP_DIR
                REPLACE_FILES )
		], @args);

  throw("need a picard_dir") if ! $picard_dir;
  $self->picard_dir($picard_dir);

  throw("need a java executable") if !$self->program;

  $self->replace_files( defined $replace_files ? $replace_files : 1);
  $self->jvm_options( defined $jvm_options ? $jvm_options : '-Xmx2g' );

  $self->tmp_dir( defined $tmp_dir ? $tmp_dir : $self->working_dir );



  $self->options($options);

  if ($index_bai && $index_bam_bai) {
      throw( "use -index_bai OR -index_bam_bai, not both" );
  }
  $self->flags('index_bam_bai', $index_bam_bai);
  $self->flags('index_bai', $index_bai);

  if (! $self->job_name) {
      $self->generate_job_name;
  }


  return $self;
}


=head2 run

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : string, either 'remove_duplicates' or 'merge'
  Function  : uses picard to process the files in $self->input_files.
  Output is files are stored in $self->output_files
  Returntype: 
  Exceptions: 
  Example   : $self->run('merge');

=cut

sub run{
    my ($self, $command) = @_;

    if ($command eq 'remove_duplicates') {
        $self->run_remove_duplicates;
    }
    elsif ($command eq 'merge') {
        $self->run_merge;
    }
    else {
        throw("Did not recognise command $command");
    }
}


=head2 run_remove_duplicates

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses MarkDuplicates.jar to remove duplicates
  Returntype: 
  Exceptions: 
  Example   : $self->run_remove_duplicates

=cut

sub run_remove_duplicates{
    my ($self) = @_;

    foreach my $input (@{$self->input_files}) {

        my $name = fileparse($input, qr/\.[sb]am/);
        my $prefix = $self->working_dir . '/' . $name . '.rmdup';
        $prefix =~ s{//}{/}g;

        my $bam = $prefix . '.bam';
        my $metrics = $prefix . '.metrics';

        my $jar = $self->picard_dir . '/MarkDuplicates.jar';
        $jar =~ s{//}{/}g;

        my $cmd = $self->program;
        $cmd .= ' ' . $self->jvm_options if ($self->jvm_options);
        $cmd .= ' -jar ' . $jar;
        $cmd .= ' ' . $self->options if ($self->options);
        $cmd .= ' TMP_DIR=' . $self->tmp_dir;
        $cmd .= ' INPUT=' . $input;
        $cmd .= ' OUTPUT=' . $bam;
        $cmd .= ' METRICS_FILE=' . $metrics;
        $cmd .= ' REMOVE_DUPLICATES=true';

        if ($self->flags('index_bai') || $self->flags('index_bam_bai')) {
            $cmd .= ' CREATE_INDEX=true';
        }

        $self->execute_command_line($cmd);

        $self->output_files($bam);
        $self->files_to_delete($metrics);

        if ($self->flags('index_bai') || $self->flags('index_bam_bai')) {
            my $bai = "$prefix.bai";
            if ($self->flags('index_bam_bai')) {
                my $bam_bai = "$prefix.bam.bai";
                $self->execute_command_line("mv $bai $bam_bai");
                $bai = $bam_bai;
            }
            $self->output_files($bai);
        }

        if ($self->replace_files) {
            $self->files_to_delete( $input );
        }
        
    }
    return;

}


=head2 run_merge

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses MergeSamFiles.jar to merge, sort and index
  Returntype: 
  Exceptions: 
  Example   : $self->run_merge

=cut


sub run_merge{
    my ($self) = @_;

    my $prefix = $self->working_dir . '/' . $self->job_name;
    $prefix =~ s{//}{/}g;

    my $bam = "$prefix.merge.bam";

    my $jar = $self->picard_dir . '/MergeSamFiles.jar';
    $jar =~ s{//}{/}g;

    my $cmd = $self->program;
    $cmd .= ' ' . $self->jvm_options if ($self->jvm_options);
    $cmd .= ' -jar ' . $jar;
    $cmd .= ' ' . $self->options if ($self->options);
    $cmd .= ' TMP_DIR=' . $self->tmp_dir;
    foreach my $input (@{$self->input_files}) {
        $cmd .= ' INPUT=' . $input;
    }
    $cmd .= ' OUTPUT=' . $bam;

    if ($self->flags('index_bai') || $self->flags('index_bam_bai')) {
        $cmd .= ' CREATE_INDEX=true';
    }

    $self->execute_command_line($cmd);

    $self->output_files($bam);

    if ($self->flags('index_bai') || $self->flags('index_bam_bai')) {
        my $bai = "$prefix.bai";
        if ($self->flags('index_bam_bai')) {
            my $bam_bai = "$prefix.bam.bai";
            $self->execute_command_line("mv $bai $bam_bai");
            $bai = $bam_bai;
        }
        $self->output_files($bai);
    }

    if ($self->replace_files) {
        $self->files_to_delete( $self->input_files );
    }


    return;
}

sub replace_files {
    my $self = shift;
    if (@_) {
        $self ->{'replace_files'} = (shift) ? 1 : 0;
    }
    return $self->{'replace_files'};
}

sub picard_dir {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'picard_dir'} = $arg;
    }
    return $self->{'picard_dir'};
}

sub jvm_options {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'jvm_options'} = $arg;
    }
    return $self->{'jvm_options'};
}

sub options {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'options'} = $arg;
    }
    return $self->{'options'};
}

sub tmp_dir {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'tmp_dir'} = $arg;
    }
    return $self->{'tmp_dir'};
}



=head2 flags

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : string, name of key e.g. "index_bai" or "index_bam_bai"
  Arg [3]   : boolean, optional, turn flag on or off
  Function  : accessor method for flags to use the various samtools utilities.
  Returntype: boolean, flag
  Exceptions: n/a
  Example   : my $run_sort_flag = $self->flags{'index_bai'};

=cut

sub flags {
    my ($self, $flag_name, $flag_value) = @_;

    if (! $self->{'flags'}) {
        $self->{'flags'} = {};
    }

    throw( "flag_name not specified")
        if (! $flag_name);

    if ($flag_value) {
        $self->{'flags'}->{$flag_name} = $flag_value ? 1 : 0;
    }

    return $self->{'flags'}->{$flag_name};
}



1;

