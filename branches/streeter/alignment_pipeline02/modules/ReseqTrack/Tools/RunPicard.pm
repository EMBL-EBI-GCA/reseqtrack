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
use ReseqTrack::Tools::FileSystemUtils qw(check_file_does_not_exist make_directory);
use File::Basename qw(fileparse);
use File::Copy qw (move);

use base qw(ReseqTrack::Tools::RunProgram);

=head2 DEFAULT_OPTIONS

  Function  : Called by the RunProgram parent object in constructor
  Returntype: hashref
  Example   : my %options = %{&ReseqTrack::Tools:RunPicard::DEFAULT_OPTIONS};

=cut

sub DEFAULT_OPTIONS { return {
        'assume_sorted' => undef, # assume input files are sorted, even if header says otherwise
        'sort_order' => 'coordinate', # can be 'coordinate' or 'queryname'
        'index_ext' => '.bam.bai', #can be '.bai' or '.bam.bai'
        'verbosity' => 'INFO',
        'max_records_in_ram' => 500000,
        'remove_duplicates' => 0, # used by run_mark_duplicates
        'use_threading' => 0, # used by merge
        'validation_stringency' => undef,
        };
}

=head2 new

  Arg [-picard_dir]   :
      string, directory containing picard jar files
  Arg [-java_exe]   :
      string, the java executable, default is 'java'
  Arg [-jvm_options]   :
      string, options for java, default is '-Xmx4g'
  Arg [-create_index]   :
      boolean, flag to create index files for all outputs
  + Arguments for ReseqTrack::Tools::RunProgram parent class

  Function  : Creates a new ReseqTrack::Tools::RunPicard object.
  Returntype: ReseqTrack::Tools::RunPicard
  Exceptions: 
  Example   : my $run_picard = ReseqTrack::Tools::RunPicard->new(
                -input_files => ['/path/sam1', '/path/sam2'],
                -picard_dir => '/path/to/picard/',
                -working_dir => '/path/to/dir/',
                -create_index => 1);

=cut


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my ( $java_exe, $picard_dir, $jvm_options, $create_index,
        )
    = rearrange( [
         qw( JAVA_EXE PICARD_DIR JVM_OPTIONS CREATE_INDEX
        )], @args);

  if (!$self->program && !$java_exe) {
    $java_exe = 'java';
  }
  $self->java_exe($java_exe);

  $self->picard_dir($picard_dir || $ENV{PICARD});
  $self->jvm_options( defined $jvm_options ? $jvm_options : '-Xmx4g' );
  $self->create_index( $create_index );

  return $self;
}


=head2 run_program

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : string, command, must be one of the following:
              'mark_duplicates', 'merge', 'sort' or 'alignment_metrics' 
  Function  : uses picard to process the files in $self->input_files.
              Output is files are stored in $self->output_files
              This method is called by the RunProgram run method
  Returntype: 
  Exceptions: Throws if it cannot find the picard_dir. Throws if command is not recognised.
  Example   : $self->run('merge');

=cut

sub run_program{
    my ($self, $command) = @_;

    throw("need a picard_dir") if ! $self->picard_dir;

    my %subs = ( 'mark_duplicates'   => \&run_mark_duplicates,
                'merge'             => \&run_merge,
                'sort'              => \&run_sort,
                'alignment_metrics' => \&run_alignment_metrics ,
    );

    throw("Did not recognise command $command") if (!defined $subs{$command});

    &{$subs{$command}}($self);
    return;
}



=head2 run_mark_duplicates

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses MarkDuplicates.jar to mark duplicates
  Returntype: 
  Exceptions: 
  Example   : $self->run_mark_duplicates

=cut

sub run_mark_duplicates {
    my $self = shift;

    my $jar = $self->_jar_path('MarkDuplicates.jar');
    foreach my $input (@{$self->input_files}) {
        my $suffix = $self->options('remove_duplicates') ? '.rmdup' : '.mrkdup';
        my $name = fileparse($input, qr/\.[sb]am/);
        my $prefix = $self->working_dir . '/' . $name . $suffix;
        $prefix =~ s{//}{/}g;
        my $bam = $prefix . '.bam';
        my $metrics = $prefix . '.metrics';

        my @cmd_words = ($self->java_exe);
        push(@cmd_words, $self->jvm_options) if ($self->jvm_options);
        push(@cmd_words, '-jar', $jar);
        push(@cmd_words, $self->_get_standard_options);
        push(@cmd_words, 'INPUT=' . $input);
        push(@cmd_words, 'OUTPUT=' . $bam);
        push(@cmd_words, 'METRICS_FILE=' . $metrics);
        push(@cmd_words, 'REMOVE_DUPLICATES='
                . ($self->options('remove_duplicates') ? 'true' : 'false'));
        push(@cmd_words, 'ASSUME_SORTED='
                . ($self->options('assume_sorted') ? 'true' : 'false')) if defined $self->options('assume_sorted');
        push(@cmd_words, 'CREATE_INDEX='
                . ($self->create_index ? 'true' : 'false'));

        my $cmd = join(' ', @cmd_words);

        $self->output_files($bam);
        $self->created_files($metrics);
        $self->created_files("$prefix.bai") if $self->create_index;

        $self->execute_command_line($cmd);

        if ($self->create_index) {
            my $bai = "$prefix.bai";
            my $index_ext = $self->options('index_ext');
            if ($index_ext && $index_ext ne '.bai') {
                my $corrected_bai = "$prefix" . $index_ext;
                move($bai, $corrected_bai)
                  or throw "move failed: $!";
                $bai = $corrected_bai;
            }
            $self->output_files($bai);
        }

    }
    return;
}

=head2 run_alignment_metrics

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Function  : uses CollectAlignmentSummaryMetrics.jar to generate alignment metrics file. Reads metrics. 
  Returntype: Collection of hashrefs. Keys described at http://picard.sourceforge.net/picard-metric-definitions.shtml#AlignmentSummaryMetrics
  Exceptions: 
  
=cut
sub run_alignment_metrics {
    my ($self) = @_;

    my $jar = $self->_jar_path('CollectAlignmentSummaryMetrics.jar');
    foreach my $input (@{$self->input_files}) {
        my ($input) = @{$self->input_files};
        my ($name,$dir) = fileparse($input, qr/\.[sb]am/);
        my $base_name = $dir.'/'. $name;
        
        my $output = $base_name.'.metrics';

        my @cmd_words = ($self->java_exe);
        push(@cmd_words, $self->jvm_options) if ($self->jvm_options);
        push(@cmd_words, '-jar', $jar);
        push(@cmd_words, 'INPUT=' . $input);
        push(@cmd_words, 'OUTPUT=' . $output);
        push(@cmd_words, 'ASSUME_SORTED='
                . ($self->options('assume_sorted') ? 'true' : 'false')) if defined $self->options('assume_sorted');
        my $cmd = join(' ', @cmd_words);

        $self->output_files($output);
        $self->execute_command_line($cmd);
    }

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

    my $input_bam_list = $self->input_files;
    throw("need more than two or more files to merge") if (@$input_bam_list <2);

    my $prefix = $self->working_dir . '/' . $self->job_name . ".merge";
    $prefix =~ s{//}{/}g;
    my $bam = "$prefix.bam";

    my $jar = $self->_jar_path('MergeSamFiles.jar');

    my $sort_order = $self->options('sort_order');
    if ($sort_order ne 'coordinate' && $sort_order ne 'queryname') {
      $sort_order = 'null';
    }

    my @cmd_words = ($self->java_exe);
    push(@cmd_words, $self->jvm_options) if ($self->jvm_options);
    push(@cmd_words, '-jar', $jar);
    push(@cmd_words, $self->_get_standard_options);
    push(@cmd_words, map {"INPUT=$_"} @$input_bam_list);
    push(@cmd_words, 'OUTPUT=' . $bam);
    push(@cmd_words, 'CREATE_INDEX='
            . ($self->create_index ? 'true' : 'false'));
    push(@cmd_words, 'SORT_ORDER=' . $sort_order);
    push(@cmd_words, 'ASSUME_SORTED='
            . ($self->options('assume_sorted') ? 'true' : 'false')) if defined $self->options('assume_sorted');
    push(@cmd_words, 'USE_THREADING='
            . ($self->options('use_threading') ? 'true' : 'false'));

    my $cmd = join(' ', @cmd_words);


    $self->output_files($bam);
    $self->created_files("$prefix.bai") if $self->create_index;

    $self->execute_command_line($cmd);

    if ($self->create_index) {
      my $bai = "$prefix.bai";
      my $index_ext = $self->options('index_ext');
      if ($index_ext && $index_ext ne '.bai') {
        my $corrected_bai = "$prefix" . $index_ext;
        move($bai, $corrected_bai)
          or throw "move failed: $!";
        $bai = $corrected_bai;
      }
      $self->output_files($bai);
    }


    return;
}

=head2 run_sort

  Arg [1]   : ReseqTrack::Tools::RunPicard
  Arg [2]   : String, sort order, defaults to coordinate
  Function  : uses SortSam.jar to sort
  Returntype: 
  Exceptions: 
  Example   : $self->run_sort

=cut

sub run_sort{
  my ($self) = @_;
    
  my $sort_order = $self->options('sort_order');
  if ($sort_order ne 'coordinate' && $sort_order ne 'queryname') {
      $sort_order = 'null';
  }
          
  my $jar = $self->_jar_path('SortSam.jar');

  foreach my $input (@{$self->input_files}) {
    my $name = fileparse($input, qr/\.[sb]am/);
    my $prefix = $self->working_dir . '/' . $name . '.sorted';
    $prefix =~ s{//}{/}g;
    my $output = $prefix . '.bam';

    my @cmd_words = ($self->java_exe);
    push(@cmd_words, $self->jvm_options) if ($self->jvm_options);
    push(@cmd_words, '-jar', $jar);
    push(@cmd_words, $self->_get_standard_options);
    push(@cmd_words, 'SORT_ORDER=' . $sort_order);
    push(@cmd_words, 'INPUT=' . $input);
    push(@cmd_words, 'OUTPUT=' . $output);
    push(@cmd_words, 'CREATE_INDEX='
            . ($self->create_index ? 'true' : 'false'));

    my $cmd = join(' ', @cmd_words);
  
    $self->output_files($output);
    $self->created_files("$prefix.bai") if $self->create_index;

    $self->execute_command_line($cmd);

    if ($self->create_index) {
        my $bai = "$prefix.bai";
        my $index_ext = $self->options('index_ext');
        if ($index_ext && $index_ext ne '.bai') {
            my $corrected_bai = "$prefix" . $index_ext;
            move($bai, $corrected_bai)
              or throw "move failed: $!";
            $bai = $corrected_bai;
        }
        $self->output_files($bai);
    }
  
  }

    return;   
}

sub _get_standard_options {
    my $self = shift;

    my @option_strings;
    push(@option_strings, 'TMP_DIR=' . $self->get_temp_dir);
    push(@option_strings, 'VERBOSITY=' . $self->options('verbosity'))
        if ($self->options('verbosity'));
    push(@option_strings, 'MAX_RECORDS_IN_RAM=' . $self->options('max_records_in_ram'))
        if ($self->options('max_records_in_ram'));
    push(@option_strings, 'VALIDATION_STRINGENCY=' . $self->options('validation_stringency'))
        if ($self->options('validation_stringency'));

    return join(' ', @option_strings);
}

sub _jar_path{
    my ($self,$jar_file) = @_;
    my $jar = $self->picard_dir . '/' . $jar_file;
    $jar =~ s{//}{/}g;

    return $jar;
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

sub java_exe {
  my $self = shift;
  return $self->program(@_);
}

sub create_index {
  my $self = shift;
  if (@_) {
    $self->{'create_index'} = (shift) ? 1 : 0;
  }
  return $self->{'create_index'};
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

