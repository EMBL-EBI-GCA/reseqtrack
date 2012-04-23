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
        $flags,
        )
    = rearrange( [
         qw( REFERENCE_INDEX REFERENCE
                FLAGS
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

  $self->flags($flags);

  return $self;
}

sub run_fix_and_calmd {
    my $self = shift;

    my $samtools = $self->program;
    my $reference = $self->reference;

    my $tmp_dir = $self->working_dir()
          .'/'.$self->job_name.'.'.$$.'.tmp/';
    check_file_does_not_exist($tmp_dir);
    $self->created_files($tmp_dir);
    make_directory($tmp_dir);

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );
      my $output_bam = $self->working_dir . "/$prefix.fixed.bam";
      $output_bam =~ s{//}{/}g;

      my @cmds;
      if ($input =~ /\.sam$/) {
        my $sam_to_bam_cmd = $self->flags('use_reference_index')
                ? "$samtools view -bSu -t " . $self->find_reference_index . " $input"
                : "$samtools view -bSu $input";
        push(@cmds, $sam_to_bam_cmd);
        push(@cmds, "$samtools sort -n -o - $tmp_dir/$prefix.nsort");
      }
      else {
        push(@cmds, "$samtools sort -n -o $input $tmp_dir/$prefix.nsort");
      }
        
      push(@cmds, "$samtools fixmate /dev/stdin /dev/stdout");
      push(@cmds, "$samtools sort -o - $tmp_dir/$prefix.csort");
      push(@cmds, "$samtools calmd -Erb - $reference");

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

    my $samtools = $self->program;

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );

      my $bam = $self->working_dir . "/$prefix.bam";
      $bam =~ s{//}{/}g;

      my $cmd = $self->program . " view -bS ";

      if ($self->flags('use_reference_index')) {
          $cmd .= "-t " . $self->find_reference_index . " ";
      }

      $cmd .= "$input > $bam";

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

    my $samtools = $self->program;

    my $tmp_dir = $self->working_dir()
          .'/'.$self->job_name.'.'.$$.'.tmp/';
    check_file_does_not_exist($tmp_dir);
    $self->created_files($tmp_dir);
    make_directory($tmp_dir);

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );
      my $output_bam = $self->working_dir . "/$prefix.sorted.bam";
      $output_bam =~ s{//}{/}g;

      my @cmds;
      if ($input =~ /\.sam$/) {
        my $sam_to_bam_cmd = $self->flags('use_reference_index')
                ? "$samtools view -bSu -t " . $self->find_reference_index . " $input"
                : "$samtools view -bSu $input";
        push(@cmds, $sam_to_bam_cmd);
        push(@cmds, "$samtools sort -n -o - $tmp_dir/$prefix.sort");
      }
      else {
        push(@cmds, "$samtools sort -n -o $input $tmp_dir/$prefix.sort");
      }

      my $cmd = join(' | ', @cmds) . " > $output_bam";

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

    my $samtools = $self->program;
    my $flag_sort = $self->flags('sort_inputs');

    my $output_bam = $self->working_dir . '/' . $self->job_name . '.merged.bam';
    $output_bam =~ s{//}{/}g;

    my $tmp_dir;
    if ($flag_sort) {
        $tmp_dir = $self->working_dir()
              .'/'.$self->job_name.'.'.$$.'.tmp/';
        check_file_does_not_exist($tmp_dir);
        $self->created_files($tmp_dir);
        make_directory($tmp_dir);
    }

    my $cmd = "$samtools merge $output_bam";

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );
      if ($input =~ /\.sam$/) {
        my $sam_to_bam_cmd = $self->flags('use_reference_index')
                ? "$samtools view -bSu -t " . $self->find_reference_index . " $input"
                : "$samtools view -bSu $input";
        if ($flag_sort) {
          $cmd .= " <($sam_to_bam_cmd | $samtools sort -n -o - $tmp_dir/$prefix.sort)";
        }
        else {
          $cmd .= " <($sam_to_bam_cmd)"
        }
      }
      else {
        if ($flag_sort) {
            $cmd .= " <($samtools sort -n -o $input $tmp_dir/$prefix.sort)";
        }
        else {
          $cmd .= " $input";
        }
      }
    }

    $self->output_files($output_bam);
    $self->execute_command_line($cmd);

}

sub run_remove_duplicates {
    my $self = shift;

    return;
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

    my %subs = {'remove_duplicates' => \&run_remove_duplicates,
                'merge'             => \&run_merge,
                'sort'              => \&run_sort,
                'index'             => \&run_index,
                'fix_and_calmd'     => \&run_fix_and_calmd ,
                'sam_to_bam'        => \&run_sam_to_bam,
    }

    throw("Did not recognise command $command") if (!defined $subs{$command});

    $self->{$subs{$command}};

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


1;

