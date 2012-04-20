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
use ReseqTrack::Tools::FileSystemUtils qw (check_file_exists);
use List::Util qw (first);
use Env qw( @PATH );
use File::Copy qw (move copy);
use Scalar::Util qw( refaddr );

use base qw(ReseqTrack::Tools::RunProgram);

=head2 new

  Arg [-reference_index]   :
      string, path of the reference genome .fai file
  Arg [-reference]   :
      string, path of the reference
  Arg [-flag_merge]   :
      boolean, flag to merge all output to a single bam
  Arg [-flag_sort]   :
      boolean, flag to sort the bam files
  Arg [-flag_index]   :
      boolean, flag to index the bam files
  Arg [-flag_sam_to_bam]   :
      boolean, flag to convert input sams to bams
  Arg [-flag_use_header]   :
      boolean, flag to use header in input sams when converting to bam (i.e. don't use an index file)
  Arg [-options_merge]   :
      string, command line options to use with "samtools merge"
  Arg [-options_sort]   :
      string, command line options to use with "samtools sort"
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
                -output_to_working_dir => 1 );

=cut


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $reference_index, $reference,
        $flag_merge, $flag_sort, $flag_index, $flag_sam_to_bam, $flag_use_header,
        $flag_calmd,
        $options_merge, $options_sort, $options_calmd,
        )
    = rearrange( [
         qw( REFERENCE_INDEX REFERENCE
                FLAG_MERGE FLAG_SORT FLAG_INDEX FLAG_SAM_TO_BAM FLAG_USE_HEADER
                FLAG_CALMD
                OPTIONS_MERGE OPTIONS_SORT OPTIONS_CALMD)
		], @args);

  #setting defaults
  if (!$self->program) {
    if ($ENV{SAMTOOLS}) {
      $self->program($ENV{SAMTOOLS} . '/samtools');
    }
    else {
      $self->program(first {-x $_} map {"$_/samtools"} @PATH);
    }
  }

  $self->reference_index($reference_index);
  $self->reference($reference);

  $self->options('merge', $options_merge);
  $self->options('sort', $options_sort);
  $self->options('calmd', $options_calmd);

  $self->flags('merge', $flag_merge);
  $self->flags('sort', $flag_sort);
  $self->flags('index', $flag_index);
  $self->flags('sam_to_bam', $flag_sam_to_bam);
  $self->flags('use_header', $flag_use_header);
  $self->flags('calmd', $flag_calmd);


  return $self;
}

sub find_reference_index {
    my ($self) = @_;

    if (! $self->reference_index) {
      my $reference_index = $self->reference;
      $reference_index =~ s/\.gz$//;
      $reference_index .= '.fai';
      $self->reference_index($reference_index);
    }

    check_file_exists($self->reference_index);
    return;
}

=head2 run_sam_to_bam

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, path of sam file
  Function  : uses samtools to convert sam to bam
  Returntype: string, path of bam file
  Exceptions: throws if there is no reference index file
  Example   : my $bam = $self->run_sam_to_bam($sam);

=cut

sub run_sam_to_bam {
    my ($self, $input_sam) = @_;

    my $prefix = fileparse($input_sam, qr/\.sam/ );

    my $bam = $self->working_dir . "/$prefix.bam";
    $bam =~ s{//}{/}g;

    my $cmd = $self->program . " view -bS ";

    if (! $self->flags('use_header')) {
        $self->find_reference_index;
        $cmd .= "-t " . $self->reference_index . " ";
    }

    $cmd .= $input_sam . " > ";
    $cmd .= $bam;

    $self->created_files($bam);
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

    my $prefix = fileparse($input_bam, qr/\.bam/ );

    my $sorted_bam_prefix = $self->working_dir . "/$prefix.srt";
    $sorted_bam_prefix =~ s{//}{/}g;

    my $sorted_bam = $sorted_bam_prefix . ".bam";

    my $cmd = $self->program . " sort";
    if ($self->options('sort')) {
        $cmd .= " " . $self->options('sort');
    }
    $cmd .= " $input_bam $sorted_bam_prefix";

    $self->created_files($sorted_bam);
    $self->execute_command_line($cmd);

    return $sorted_bam;
}

sub run_calmd {
    my ($self, $input_bam) = @_;

    my $prefix = fileparse($input_bam, qr/\.bam/ );

    my $calmd_bam = $self->working_dir . "/$prefix.calmd.bam";
    $calmd_bam =~ s{//}{/}g;

    my $cmd = $self->program . " calmd";

    my $options = defined $self->options('calmd') ? $self->options('calmd') : '-Erb';
    $cmd .= " " . $options;

    $cmd .= " $input_bam ";
    $cmd .= $self->reference;
    $cmd .= " > $calmd_bam";

    $self->created_files($calmd_bam);
    $self->execute_command_line($cmd);

    return $calmd_bam;
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

    $self->created_files($output_bai);
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

    my $merged_bam = $self->working_dir . '/' . $self->job_name . '.merged.bam';
    $merged_bam =~ s{//}{/}g;
    $self->created_files($merged_bam);

    if (@$input_bam_list >1) {
      my $cmd = $self->program . " merge";
      if ($self->options('merge')) {
          $cmd .= " " . $self->options('merge');
      }
      $cmd .= " " . $merged_bam;
      foreach my $bam (@$input_bam_list) {
          $cmd .= " " . $bam;
      }
      $self->execute_command_line($cmd);
    }
    elsif (@$input_bam_list ==1) {
      my $input_bam = $input_bam_list->[0];
      if (grep {$_ eq $input_bam} @{$self->created_files}) {
        print "renaming $input_bam to $merged_bam\n";
        move($input_bam, $merged_bam)
              or throw "copy failed: $!";
      }
      else {
        print "copying $input_bam to $merged_bam\n";
        copy($input_bam, $merged_bam)
              or throw "copy failed: $!";
      }
    }

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

sub run_program {
    my ($self) = @_;

    my $input_files = $self->input_files;
    my $current_files = $input_files;

    if ($self->flags('sam_to_bam')) {
        my @bams;
        foreach my $sam (@$current_files) {
            my $bam = $self->run_sam_to_bam($sam);
            push(@bams, $bam);
        }

        $current_files = \@bams;
    }

    if ($self->flags('sort')) {
        my @sorted_bams;
        foreach my $file (@$current_files) {
            my $sorted_bam = $self->run_sort($file);
            push(@sorted_bams, $sorted_bam);
        }
        
        $current_files = \@sorted_bams;
    }

    if ($self->flags('merge')) {
        my $merged_bam = $self->run_merge($current_files);
        $current_files = [$merged_bam];
    }

    if ($self->flags('calmd')) {
      my @calmd_bams;
      foreach my $file(@$current_files) {
        my $calmd_bam = $self->run_calmd($file);
        push(@calmd_bams, $calmd_bam);
      }
      $current_files = \@calmd_bams;
    }

    if (refaddr($input_files) != refaddr($current_files)) {
      $self->output_files($current_files);
    }

    if ($self->flags('index')) {
        foreach my $file (@$current_files) {
            my $index = $self->run_index($file);
            $self->output_files($index);
        }
    }

    return;
}

=head2 reference

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, optional, path of genome reference
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

=head2 reference_index

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, optional, path of reference fai file
  Function  : accessor method for reference_index
  Returntype: string
  Exceptions: n/a
  Example   : my $reference_index = $self->reference_index;

=cut

sub reference_index {
  my ($self, $reference_index) = @_;
  if ($reference_index) {
    $self->{'reference_index'} = $reference_index;
  }
  return $self->{'reference_index'};
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

sub output_bai_files {
  my $self = shift;
  my @files = grep {/\.bai$/} @{$self->output_files};
  return \@files;
}

sub output_bam_files {
  my $self = shift;
  my @files = grep {/\.bam$/} @{$self->output_files};
  return \@files;
}

1;

