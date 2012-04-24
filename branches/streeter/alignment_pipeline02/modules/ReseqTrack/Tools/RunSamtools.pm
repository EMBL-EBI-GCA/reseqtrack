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
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_file_does_not_exist make_directory);
use List::Util qw (first);
use Env qw( @PATH );

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

sub DEFAULT_OPTIONS { return [
        'use_reference_index => 0,
        'input_sort_status' => '', # can be 'c' for coordinate or 'n' for name
        'force_overwrite' => 1,
        'attach_RG_tag' => 0,
        'compute_BQ_tag' => 1,
        'ext_BAQ_calc' => 1,
        ];
}

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $reference_index, $reference,
        )
    = rearrange( [
         qw( REFERENCE_INDEX REFERENCE
                ) ], @args);

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

  return $self;
}

sub run_fix_and_calmd {
    my $self = shift;

    foreach my $input (@{$self->input_files}) {
      my $output_bam = $self->working_dir . "/$prefix.fixed.bam";
      $output_bam =~ s{//}{/}g;

      my @cmds;
      push(@cmds, $self->_get_file_to_sorted_bam_cmd($input, 1, 1);
      push(@cmds, $self->_get_fixmate_cmd);
      push(@cmds, $self->_get_sort_cmd);
      push(@cmds, $self->_get_calmd_cmd(1));

      my $cmd = join(' | ', @cmds) . " > $output_bam";

      $self->output_files($output_bam);
      $self->execute_command_line($cmd);
    }
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
    return $self->reference_index;
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
    my ($self) = @_;

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );

      my $bam = $self->working_dir . "/$prefix.bam";
      $bam =~ s{//}{/}g;

      my $cmd = $self->_get_sam_to_bam_cmd($input);
      $cmd .= " > $bam";

      $self->created_files($bam);
      $self->execute_command_line($cmd);
    }

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
    my $self = shift;

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );
      my $output_bam = $self->working_dir . "/$prefix.sorted.bam";
      $output_bam =~ s{//}{/}g;

      my $cmd = $self->_get_file_to_sorted_bam_cmd($input);
      $cmd .= " > $output_bam";

      $self->output_files($output_bam);
      $self->execute_command_line($cmd);
    }
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
    my $self = shift;

    foreach my $file (@{$self->input_files}) {
        my $output_bai = $file . ".bai";
        my $cmd = $self->program . " index " . $file;

        $self->output_files($output_bai);
        $self->execute_command_line($cmd);
    }
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
    my $self = shift;

    throw("need more than two or more files to merge") if (@{$self->input_files} <2);

    my $output_bam = $self->working_dir . '/' . $self->job_name . '.merged.bam';
    $output_bam =~ s{//}{/}g;
    throw("file already exists and force_overwrite is not set")
            if (-e $output_bam && ! $self->options('force_overwrite');

    my @cmd_words = ($self->program, 'merge');
    push(@cmd_words, '-f') if ($self->options('force_overwrite'));
    push(@cmd_words, '-r') if ($self->options('attach_RG_tag'));
    push(@cmd_words, $output_bam);

    foreach my $input (@{$self->input_files}) {
      if ($self->options('input_sort_status') ne 'c' || $input =~ /\.sam$/) {
        my $file_to_sorted_bam_cmd = $self->_get_file_to_sorted_bam_cmd($input, 0, 1);
        push(@cmd_words, '<(', $sam_to_sorted_bam_cmd, ')');
      }
      else {
          push(@cmd_words, $input);
      }
    }
    my $cmd = (join(' '), @cmd_words);

    $self->output_files($output_bam);
    $self->execute_command_line($cmd);

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
    my ($self, $command) = @_;

    my %subs = {'merge'             => \&run_merge,
                'sort'              => \&run_sort,
                'index'             => \&run_index,
                'fix_and_calmd'     => \&run_fix_and_calmd ,
                'sam_to_bam'        => \&run_sam_to_bam,
    };

    throw("Did not recognise command $command") if (!defined $subs{$command});

    &{$subs{$command}}($self);

    return;
}



sub _get_fixmate_cmd {
    my ($self, $bam) = @_;
    my @cmd_words = ($self->program, 'fixmate');
    push(@cmd_words, $bam || '/dev/stdin');
    push(@cmd_words, '/dev/stdout');
    return join(' ', @cmd_words);
}

sub _get_calmd_cmd {
    my ($self, $uncompressed, $input) = @_;
    my @cmd_words = ($self->program, 'calmd');
    push(@cmd_words, $uncompressed ? '-u' : '-b');
    push(@cmd_words, '-r') if ($self->options('compute_BQ_tag');
    push(@cmd_words, '-E') if ($self->options('ext_BAQ_calc');

    if ($input) {
      push(@cmd_words, '-S') if ($input =~ /\.sam$/);
      push(@cmd_words, $input);
    }
    else {
      push(@cmd_words, '-');
    }
    push(@cmd_words, $self->reference);
    return join(' ', @cmd_words);
}

sub _get_file_to_sorted_bam_cmd {
    my ($self, $file, $name_sort, $uncompressed) @_;

    return $self->_get_sam_to_sorted_bam_cmd($file, $name_sort, $uncompressed)
        if ($input =~ /\.sam$/);

    return $self->_get_sort_cmd(1, $file)
        if ($name_sort && $self->options('input_sort_status') ne 'n');

    return $self->_get_sort_cmd(0, $file)
        if ($self->options('input_sort_status') ne 'c');

    my @cmd_words = ($self->program, 'view', '-hb');
    push(@cmd_words, '-u') if $uncompressed;
    push(@cmd_words, $file);
    return;
}


sub _get_sam_to_sorted_bam_cmd {
    my ($self, $sam, $name_sort, $uncompressed) = @_;

    my $pipe_to_sort = $name_sort && $self->options('input_sort_status) ne 'n';
    $pipe_to_sort ||= !$name_sort && $self->options('input_sort_status) ne 'c';

    my @cmds = ($self->_get_sam_to_bam_cmd($sam, $uncompressed || $pipe_to_sort));
    if ($pipe_to_sort) {
        push(@cmds, $self->_get_sort_cmd($name_sort));
    }
    return join(' | ', @cmds);
}

sub _get_sam_to_bam_cmd {
    my ($self, $sam, $uncompressed) = @_;

    my @cmd_words = ($self->program, 'view', '-bS');
    push(@cmd_words, '-u') if ($uncompressed);
    push(@cmd_words, '-t', $self->find_reference_index)
            if ($self->options('use_reference_index'));
    push(@cmd_words, $sam || '-');

    return join(' ', @cmd_words);
}

sub _get_sort_cmd {
    my ($self, $name_sort, $bam) = @_;

    my $prefix = $self->get_temp_dir .'/'. $self->job_name;
    $prefix .= '.' . $name_sort ? 'nsort' : 'csort';

    my @cmd_words = ($self->program, 'sort', '-o');
    push(@cmd_words, '-n') if ($name_sort);
    push(@cmd_words, $bam || '-');
    push(@cmd_words, $prefix);

    return join(' ', @cmd_words);
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



1;

