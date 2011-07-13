=pod

=head1 NAME

ReseqTrack::Tools::RunSamtools

=head1 SYNOPSIS

This is a class for running samtools
It is a sub class of a ReseqTrack::Tools::RunProgram.

=head1 Example


=cut

package ReseqTrack::Tools::RunSamtools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-flag_merge]   :
      boolean, flag to merge all output to a single bam
  Arg [-flag_sort]   :
      boolean, flag to sort the bam files
  Arg [-flag_index]   :
      boolean, flag to index the bam files
  Arg [-flag_sam_to_bam]   :
      boolean, flag to convert input sams to bams
  Arg [-options_merge]   :
      string, command line options to use with "samtools merge"
  Arg [-options_sort]   :
      string, command line options to use with "samtools sort"
  Arg [-replace_files]   :
      boolean, default 1, any input / intermediate file will be deleted after output is created.
  Arg [-output_to_working_dir]   :
      boolean, default 0, flag to write output files to the working directory rather than the input directory
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunSamtools object.
  Returntype: ReseqTrack::Tools::RunSamtools
  Exceptions: 
  Example   : my $run_alignment = ReseqTrack::Tools::RunSamtools->new(
                -input_files => ['/path/sam1', '/path/sam2'],
                -program => "samtools",
                -working_dir => '/path/to/dir/',
                -flag_sam_to_bam => 1,
                -flag_merge => 1,
                -flag_sort => 1,
                -flag_index => 1,
                -options_merge => "-r",
                -options_sort => "-m 100000000"
                -replace_files => 1,
                -output_to_working_dir => 1 );

=cut


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $reference,
        $flag_merge, $flag_sort, $flag_index, $flag_sam_to_bam,
        $options_merge, $options_sort,
        $replace_files, $output_to_working_dir)
    = rearrange( [
         qw( REFERENCE
                FLAG_MERGE FLAG_SORT FLAG_INDEX FLAG_SAM_TO_BAM
                OPTIONS_MERGE OPTIONS_SORT
                REPLACE_FILES OUTPUT_TO_WORKING_DIR )
		], @args);

  $self->reference($reference);
  $self->output_to_working_dir($output_to_working_dir);

  if (! defined $replace_files) {
      $replace_files = 1;
  }
  $self->replace_files($replace_files);

  $self->options('merge', $options_merge);
  $self->options('sort', $options_sort);

  $self->flags('merge', $flag_merge);
  $self->flags('sort', $flag_sort);
  $self->flags('index', $flag_index);
  $self->flags('sam_to_bam', $flag_sam_to_bam);

  if (! $self->job_name) {
      $self->generate_job_name;
  }


  return $self;
}

=head2 run_sam_to_bam

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, path of sam file
  Function  : uses samtools to convert sam to bam
  Returntype: string, path of bam file
  Exceptions: 
  Example   : my $bam = $self->run_sam_to_bam($sam);
=cut

sub run_sam_to_bam {
    my ($self, $input_sam) = @_;

    my ($prefix, $dir) = fileparse($input_sam, qr/\.sam/ );

    if ($self->output_to_working_dir) {
        $dir = $self->working_dir;
    }

    my $bam = "$dir/$prefix.bam";
    $bam =~ s{//}{/}g;

    my $cmd = $self->program . " import ";
    $cmd .= $self->reference . " ";
    $cmd .= $input_sam . " ";
    $cmd .= $bam;

    $self->execute_command_line($cmd);

    return $bam;
}

=head2 run_sort

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, path of bam file
  Function  : uses samtools to sort bam file
  Returntype: string, path of sorted bam file
  Exceptions: 
  Example   : my $sorted_bam = $self->run_sort($bam);
=cut

sub run_sort {
    my ($self, $input_bam) = @_;

    my ($prefix, $dir) = fileparse($input_bam, qr/\.bam/ );

    if ($self->output_to_working_dir) {
        $dir = $self->working_dir;
    }

    my $sorted_bam_prefix = "$dir/$prefix.srt";
    $sorted_bam_prefix =~ s{//}{/}g;

    my $sorted_bam = $sorted_bam_prefix . ".bam";

    my $cmd = $self->program . " sort";
    if ($self->options('sort')) {
        $cmd .= " " . $self->options('sort');
    }
    $cmd .= " $input_bam $sorted_bam_prefix";

    $self->execute_command_line($cmd);

    return $sorted_bam;
}


=head2 run_index

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, path of bam file
  Function  : uses samtools to index bam file
  Returntype: string, path of index file
  Exceptions: 
  Example   : my $bai = $self->run_index($bam);
=cut

sub run_index {
    my ($self, $input_bam) = @_;

    my $output_bai = $input_bam . ".bai";

    my $cmd = $self->program . " index " . $input_bam;

    $self->execute_command_line($cmd);

    return $output_bai;
}


=head2 run_merge

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : arrayref of strings, paths of bam files
  Function  : uses samtools to merge bam files
  Returntype: string, path of merged bam file
  Exceptions: 
  Example   : my $merged_bam = $self->run_index([$bam1, $bam2]);
=cut

sub run_merge {
    my ($self, $input_bam_list) = @_;

    my $dir = $self->output_to_working_dir
                            ? $self->working_dir
                            : [fileparse( $input_bam_list->[0] )]->[1];

    my $merged_bam = $dir . '/' . $self->job_name . '_merged.bam';
    $merged_bam =~ s{//}{/}g;

    my $cmd = $self->program . " merge";
    if ($self->options('merge')) {
        $cmd .= " " . $self->options('merge');
    }
    $cmd .= " " . $merged_bam;
    foreach my $bam (@$input_bam_list) {
        $cmd .= " " . $bam;
    }

    $self->execute_command_line($cmd);

    return $merged_bam;
}



=head2 run

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Function  : uses samtools to process the files in $self->input_files.
  Output is files are stored in $self->output_files
  Returntype: 
  Exceptions: 
  Example   : $self->run();
=cut

sub run {
    my ($self) = @_;

    my @current_files = @{$self->input_files};
    my @intermediate_files;
    my $current_files_are_input = 1;

    if ($self->flags('sam_to_bam')) {
        my @bams;
        foreach my $sam (@current_files) {
            my $bam = $self->run_sam_to_bam($sam);
            push(@bams, $bam);
        }

        @current_files = @bams;
        $current_files_are_input = 0;
    }

    if ($self->flags('sort')) {
        my @sorted_bams;
        foreach my $file (@current_files) {
            my $sorted_bam = $self->run_sort($file);
            push(@sorted_bams, $sorted_bam);
        }
        
        if (! $current_files_are_input) {
            push(@intermediate_files, @current_files);
        }

        @current_files = @sorted_bams;
        $current_files_are_input = 0;
    }

    if ($self->flags('merge') && @current_files > 1) {
        my $merged_bam = $self->run_merge(\@current_files);

        if (! $current_files_are_input) {
            push(@intermediate_files, @current_files);
        }
        
        @current_files = ($merged_bam);
        $current_files_are_input = 0;
    }

    if (! $current_files_are_input) {
        if ($self->replace_files) {
            $self->files_to_delete( $self->input_files );
            $self->files_to_delete( \@intermediate_files );
            $self->output_files(\@current_files);
        }
        else {
            $self->output_files( \@intermediate_files );
            $self->output_files(\@current_files);
        }
    }


    if ($self->flags('index')) {
        foreach my $file (@current_files) {
            my $index = $self->run_index($file);
            $self->output_files($index);
        }
    }

    return;

}

=head2 reference

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, optional, path of reference genome
  Function  : accessor method for reference
  Returntype: string
  Exceptions: n/a
  Example   : my $reference = $self->reference;

=cut

sub reference {
  my ($self, $reference) = @_;
  if ($reference) {
    $self->{'reference'} = $reference;
  }
  return $self->{'reference'};
}


=head2 replace_files

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : boolean, optional, value of replace_files flag
  Function  : accessor method for replace_files flag
  Returntype: boolean, replace_files flag
  Exceptions: n/a
  Example   : $self->replace_files(1);

=cut

sub replace_files {
    my $self = shift;
    
    if (@_) {
        $self ->{'replace_files'} = (shift) ? 1 : 0;
    }
    return $self->{'replace_files'};
}

=head2 output_to_working_dir

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : boolean, optional, value of output_to_working_dir flag
  Function  : accessor method for output_to_working_dir flag
  Returntype: boolean, output_to_working_dir flag
  Exceptions: n/a
  Example   : $self->output_to_working_dir(1);

=cut

sub output_to_working_dir {
    my $self = shift;
    
    if (@_) {
        $self ->{'output_to_working_dir'} = (shift) ? 1 : 0;
    }
    return $self->{'output_to_working_dir'};
}

=head2 options

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, name of key e.g. "merge" or "sort"
  Arg [3]   : string, optional, options to be used on command line e.g. "-m 100000000"
  Function  : accessor method for command line options
  Returntype: string, command line options
  Exceptions: n/a
  Example   : my $sort_options = $self->options{'sort'};

=cut

sub options {
    my ($self, $option_name, $option_value) = @_;

    if (! $self->{'options'}) {
        $self->{'options'} = {};
    }

    throw( "option_name not specified")
        if (! $option_name);

    if ($option_value) {
        $self->{'options'}->{$option_name} = $option_value ? 1 : 0;
    }

    return $self->{'options'}->{$option_name};
}

=head2 flags

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, name of key e.g. "merge" or "sort"
  Arg [3]   : boolean, optional, turn flag on or off
  Function  : accessor method for flags to use the various samtools utilities.
  Returntype: boolean, flag
  Exceptions: n/a
  Example   : my $run_sort_flag = $self->flags{'sort'};

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

