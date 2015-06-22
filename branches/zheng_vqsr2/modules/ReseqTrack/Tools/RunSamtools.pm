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
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists make_directory);

use base qw(ReseqTrack::Tools::RunProgram);


=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunSamtools::DEFAULT_OPTIONS};

=cut

sub DEFAULT_OPTIONS { return {
        'use_reference_index' => 0, # use the -t option when importing a sam
        'input_sort_status' => undef, # can be 'c' for coordinate or 'n' for name
        'output_sort_status' => undef, # can be 'c' for coordinate or 'n' for name
        'uncompressed_output' => 0,
        'force_overwrite' => 1, # used by samtools merge
        'attach_RG_tag' => 0, # used by samtools merge
        'compute_BQ_tag' => 1, # used by samtools calmd
        'ext_BAQ_calc' => 1, # used by samtools calmd
        };
}

=head2 new

  Arg [-reference_index]   :
      string, path of the reference genome .fai file
  Arg [-reference]   :
      string, path of the reference
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunSamtools object.
  Returntype: ReseqTrack::Tools::RunSamtools
  Exceptions: 
  Example   : my $run_alignment = ReseqTrack::Tools::RunSamtools->new(
                -input_files => ['/path/sam1', '/path/sam2'],
                -program => "samtools",
                -working_dir => '/path/to/dir/');

=cut

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
      $self->program('samtools');
    }
  }

  $self->reference_index($reference_index);
  $self->reference($reference);

  return $self;
}

sub run_fix_and_calmd {
    my $self = shift;

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );
      my $output_bam = $self->working_dir . "/$prefix.fixed.md.bam";
      $output_bam =~ s{//}{/}g;

      my @cmds;
      push(@cmds, $self->_get_file_to_sorted_bam_cmd($input, 1, 1));
      push(@cmds, $self->_get_fixmate_cmd);
      push(@cmds, $self->_get_sort_cmd);
      push(@cmds, $self->_get_calmd_cmd($self->options('uncompessed_output')));

      my $cmd = join(' | ', @cmds) . " > $output_bam";

      $self->output_files($output_bam);
      $self->execute_command_line($cmd);
    }
}

sub run_calmd {
    my $self = shift;

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );
      my $output_bam = $self->working_dir . "/$prefix.md.bam";
      $output_bam =~ s{//}{/}g;

      my @cmds;
      push(@cmds, $self->_get_file_to_sorted_bam_cmd($input, 0, 1));
      my $output_sort_status = $self->options('output_sort_status');
      if ($output_sort_status && $output_sort_status eq 'n') {
        push(@cmds, $self->_get_calmd_cmd(1));
        push(@cmds, $self->_get_sort_cmd(1));
      }
      else {
        push(@cmds, $self->_get_calmd_cmd($self->options('uncompressed_output')));
      }

      my $cmd = join(' | ', @cmds) . " > $output_bam";

      $self->output_files($output_bam);
      $self->execute_command_line($cmd);
    }
}

sub run_fixmate {
    my $self = shift;

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );
      my $output_bam = $self->working_dir . "/$prefix.fixed.bam";
      $output_bam =~ s{//}{/}g;

      my @cmds;
      push(@cmds, $self->_get_file_to_sorted_bam_cmd($input, 1, 1));
      push(@cmds, $self->_get_fixmate_cmd);
      if ($self->options('output_sort_status') eq 'c') {
        push(@cmds, $self->_get_sort_cmd());
      }

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

sub run_sam_to_bam {
    my ($self) = @_;

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );

      my $bam = $self->working_dir . "/$prefix.bam";
      $bam =~ s{//}{/}g;

      my $output_sort_status = $self->options('output_sort_status');

      my $cmd;
      if ($output_sort_status eq 'c') {
        $cmd = $self->_get_file_to_sorted_bam_cmd($input, 0, $self->options('uncompressed_output'));
      }
      elsif ($output_sort_status eq 'n') {
        $cmd = $self->_get_file_to_sorted_bam_cmd($input, 1, $self->options('uncompressed_output'));
      }
      else {
        $cmd = $self->_get_sam_to_bam_cmd($input, $self->options('uncompressed_output'));
      }
      $cmd .= " > $bam";

      $self->created_files($bam);
      $self->execute_command_line($cmd);
    }

}

sub run_sort {
    my $self = shift;

    foreach my $input (@{$self->input_files}) {
      my $prefix = fileparse($input, qr/\.[sb]am/ );
      my $output_bam = $self->working_dir . "/$prefix.sorted.bam";
      $output_bam =~ s{//}{/}g;

      my $cmd = $self->options('output_sort_status') eq 'n'
            ? $self->_get_file_to_sorted_bam_cmd($input, 1)
            : $self->_get_file_to_sorted_bam_cmd($input);
      $cmd .= " > $output_bam";

      $self->output_files($output_bam);
      $self->execute_command_line($cmd);
    }
}


sub run_index {
    my $self = shift;

    foreach my $file (@{$self->input_files}) {
        my $output_bai = $file . ".bai";
        my $cmd = $self->program . " index " . $file;

        $self->output_files($output_bai);
        $self->execute_command_line($cmd);
    }
}

sub run_merge {
    my $self = shift;

    throw("need more than two or more files to merge") if (@{$self->input_files} <2);

    my $output_bam = $self->working_dir . '/' . $self->job_name . '.merged.bam';
    $output_bam =~ s{//}{/}g;
    throw("file already exists and force_overwrite is not set")
            if (-e $output_bam && ! $self->options('force_overwrite'));

    my $output_sort_status = $self->options('output_sort_status');
    my $input_sort_status = $self->options('input_sort_status');

    my @cmd_words = ("bash -c '");
    push(@cmd_words, $self->program, 'merge');
    push(@cmd_words, '-f') if ($self->options('force_overwrite'));
    push(@cmd_words, '-r') if ($self->options('attach_RG_tag'));
    push(@cmd_words, '-n') if ($output_sort_status eq 'n');

    push(@cmd_words, $output_bam);

    my $needs_sorting = ($output_sort_status eq 'n' && $input_sort_status ne 'n');
    $needs_sorting ||= ($output_sort_status ne 'n' && $input_sort_status ne 'c');

    foreach my $input (@{$self->input_files}) {
      if ($needs_sorting || $input =~ /\.sam$/) {
        my $file_to_sorted_bam_cmd = $self->_get_file_to_sorted_bam_cmd($input, ($output_sort_status eq 'n'), 1);
        push(@cmd_words, '<(', $file_to_sorted_bam_cmd, ')');
      }
      else {
          push(@cmd_words, $input);
      }
    }
    push(@cmd_words, "'");
    my $cmd = join(' ', @cmd_words);

    $self->output_files($output_bam);
    $self->execute_command_line($cmd);

}


=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunSamtools
  Arg [2]   : string, command, must be one of the following:
              merge, sort, index, fix_and_calmd, calmd, sam_to_bam
  Function  : uses samtools to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if the command is not recognised.
  Example   : $self->run();

=cut

sub run_program {
    my ($self, $command) = @_;

    my %subs = ('merge'             => \&run_merge,
                'sort'              => \&run_sort,
                'index'             => \&run_index,
                'fix_and_calmd'     => \&run_fix_and_calmd ,
                'calmd'             => \&run_calmd ,
                'sam_to_bam'        => \&run_sam_to_bam,
    );

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
    push(@cmd_words, '-r') if ($self->options('compute_BQ_tag'));
    push(@cmd_words, '-E') if ($self->options('ext_BAQ_calc'));

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
    my ($self, $file, $name_sort, $uncompressed) = @_;

    return $self->_get_sam_to_sorted_bam_cmd($file, $name_sort, $uncompressed)
        if ($file =~ /\.sam$/);

    my $input_sort_status = $self->options('input_sort_status');
    return $self->_get_sort_cmd(1, $file)
        if ($name_sort && (!$input_sort_status || $input_sort_status ne 'n'));

    return $self->_get_sort_cmd(0, $file)
        if (!$input_sort_status || $input_sort_status ne 'c');

    my @cmd_words = ($self->program, 'view', '-hb');
    push(@cmd_words, '-u') if $uncompressed;
    push(@cmd_words, $file);
    return join(' ', @cmd_words);
}


sub _get_sam_to_sorted_bam_cmd {
    my ($self, $sam, $name_sort, $uncompressed) = @_;

    my $pipe_to_sort = $name_sort && $self->options('input_sort_status') ne 'n';
    $pipe_to_sort ||= !$name_sort && $self->options('input_sort_status') ne 'c';

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
    $prefix .= '.' . ($name_sort ? 'nsort' : 'csort');

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

