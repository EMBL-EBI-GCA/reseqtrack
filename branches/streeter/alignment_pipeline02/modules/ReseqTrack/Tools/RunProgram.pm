=head1 NAME

ReseqTrack::Tools::RunProgram

=head1 SYNOPSIS

This is a base class for RunProgram objects and provides some standard
accessor methods and throws exceptions when vital methods aren't implemented in
the child classes. The Child classes should wrap specific algorithms.

=head1 Example

=cut

package ReseqTrack::Tools::RunProgram;

use strict;
use warnings;

use File::Path qw(mkpath);
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::GeneralUtils qw(execute_system_command);


=head2 new

  Arg [-input_files]   :
      arrayref of strings,  paths to any input files
  Arg [-program]   :
      string, the program (if in $PATH) or path to the program
  Arg [-job_name]   :
      string, name of the job
  Arg [-working_dir]   :
      string, default '/tmp/', path of the working directory

  Arg [-echo_cmd_line]   :
      boolean flag, default 0, echo the command line rather than execute commands
  Arg [-save_files_for_deletion]   :
      boolean flag, default 0, save temporary files (i.e. do not delete them)

  Function  : Creates a new ReseqTrack::Tools::RunProgram object 
  Returntype: ReseqTrack::Tools::RunProgram
  Exceptions: 
  Example   : my $run_program = ReseqTrack::Tools::RunProgram->new(
                -input_files => ['/path/file1', '/path/file2'],
                -program => "program",
                -working_dir => '/path/to/dir/',
                -job_name => "my_job",
                -echo_cmd_line => 0,
                -save_files_for_deletion => 0 );

=cut


# package variables:
my $GLOBAL_ECHO_CMD_LINE;
my $GLOBAL_SAVE_FILES_FOR_DELETION;
my $GLOBAL_WORKING_DIR = "/tmp/";

sub new {
  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;

  #input files, we will allow single file name or list of file names
  my ( $input_files,  $program, $job_name, $working_dir,
        $options,
        $echo_cmd_line, $save_files_for_deletion )
    = rearrange( [
         qw( INPUT_FILES PROGRAM JOB_NAME WORKING_DIR
             OPTIONS
             ECHO_CMD_LINE SAVE_FILES_FOR_DELETION )
		], @args);

  $self->input_files($input_files);
  $self->program($program);
  $self->job_name($job_name);
  $self->working_dir($working_dir);
  $self->echo_cmd_line($echo_cmd_line);
  $self->save_files_for_deletion($save_files_for_deletion);

  while (my ($option_name, $option_value) = each %{$self->DEFAULT_OPTIONS}) {
      $options{$option_name} = $option_value if (! defined $options{$option_name});
  }
  while (my ($option_name, $option_value) = each %$options) {
      $self->options($option_name, $option_value);
  }

  return $self;
}

=head2 run

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Function  : each child object should implement a run method
  Returntype: 
  Exceptions: throws as this method should be implemented in the child class
  Example   : 

=cut

sub run {
  my ($self, @args) = @_;

  $self->_running(1);

  throw "do not have a program executable\n" if (! $self->program);
  check_file_exists($self->program);
  foreach my $file (@{$self->input_files}) {
    check_file_exists($file);
  }
  if (! $self->job_name) {
      $self->generate_job_name;
  }

  $self->change_dir;
  $self->run_program(@args);

  $self->_running(0);
  $self->delete_files;
}

sub DESTROY {
  my $self = shift;
  $self->delete_files;
}

=head2 execute_command_line

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, command line to execute
  Function  : executes command
  Returntype: exit code
  Exceptions: throws if execution of command fails
  Example   : $self->execute_command_line('/path/to/program -options > output');

=cut

sub execute_command_line {
    my ($self, $command_line) = @_;
    my $exit;
    if ($self->echo_cmd_line) {
        $command_line = "echo \'" . $command_line . "\'";
    }
    print $command_line . "\n";
    $exit = execute_system_command( $command_line );
    return $exit;
}

=head2 save_files_for_deletion

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : boolean, optional, value of save_files_for_deletion flag
  Function  : accessor method for save_files_for_deletion flag.  Files marked for deletion
  will not be deleted.
  Returntype: boolean
  Exceptions: 
  Example   : my $flag = $self->save_files_for_deletion;

=cut

sub save_files_for_deletion {
  my $self = shift;
  if (@_) {
    $GLOBAL_SAVE_FILES_FOR_DELETION = (shift) ? 1 : 0;
  }
  return $GLOBAL_SAVE_FILES_FOR_DELETION;
}

=head2 echo_cmd_line

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : boolean, optional, value of echo_cmd_line flag
  Function  : accessor method for echo_cmd_line flag.  Command lines will be 'echoed'
  but not executed.
  Returntype: boolean
  Exceptions: 
  Example   : my $flag = $self->echo_cmd_line;

=cut

sub echo_cmd_line {
  my $self = shift;
  if (@_) {
    $GLOBAL_ECHO_CMD_LINE = (shift) ? 1 : 0;
  }
  return $GLOBAL_ECHO_CMD_LINE;
}

=head2 program

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, optional, program (if in $PATH) or path to program
  Function  : accessor method for the program
  Returntype: string
  Exceptions: throws if executable does not exist
  Example   : my $program = $self->program;

=cut

sub program {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'program'} = $arg;
  }
  return $self->{'program'};
}

=head2 job_name

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, optional, job name
  Function  : accessor method for the job name
  Returntype: string
  Exceptions: 
  Example   : my $job_name = $self->job_name;

=cut

sub job_name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'job_name'} = $arg;
  }
  return $self->{'job_name'};
}


=head2 generate_job_name

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Function  : sets the job name equal to the process number ($$).  Child classes may
  override this function.
  Returntype: string
  Exceptions: 
  Example   : $self->generate_job_name();

=cut

sub generate_job_name {
    my $self = shift;
    return $self->job_name($$);
}


=head2 working_dir

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, path of working directory
  Function  : accessor method for working directory
  Returntype: string
  Exceptions: n/a
  Example   : my $work_dir = $self->working_dir;

=cut

sub working_dir {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $GLOBAL_WORKING_DIR = $arg;
  }
  return $GLOBAL_WORKING_DIR;
}

=head2 change_dir

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, path of working directory.  Default is $self->working_dir
  Function  : changes to given working directory.  Also creates directory
  if it does not already exist.
  Returntype: string
  Exceptions: throws if it can't change to given directory
  Example   : $self->change_dir();

=cut

sub change_dir {
  my ( $self, $dir ) = @_;

  if (! $dir){
    $dir = $self->working_dir;
  }
  check_directory_exists($dir);

  chdir($dir)
    or throw( "Failed to change to $dir");
  return $dir;
}

=head2 output_files

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for output files
  Returntype: arrayref of strings
  Exceptions: n/a
  Example   : $self->output_files('path/to/file');

=cut

sub output_files {
  my ( $self, $arg ) = @_;

  $self->{'output_files'} ||= {};
  if ($arg) {
    foreach my $file (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      $file =~ s{//}{/}g;
      $self->{'output_files'}->{$file} = 1;
    }
    $self->created_files($arg, 1);
  }

  my @files = keys %{$self->{'output_files'}};
  return \@files;
}

=head2 input_files

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for input files.  Also checks that each file exists.
  Returntype: arrayref of strings
  Exceptions: throws if file does not exist
  Example   : $self->input_files('path/to/file');

=cut

sub input_files {
  my ( $self, $arg ) = @_;

  $self->{'input_files'} ||= {};
  if ($arg) {
    foreach my $file (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      $file =~ s{//}{/}g;
      $self->{'input_files'}->{$file} = 1;
    }
  }

  my @files = keys %{$self->{'input_files'}};
  return \@files;
}


=head2 files_to_delete

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string or arrayref of strings
  Function  : store files and directories to be deleted when delete_files is called
  or when object goes out of scope
  Returntype: arrayref of strings
  Exceptions: n/a
  Example   : $self->files_to_delete('path/to/file');

=cut

sub created_files {
  my ( $self, $arg, $save_from_deletion ) = @_;

  $self->{'created_files'} ||= {};
  if ($arg) {
    $save_from_deletion = $save_from_deletion ? 1 : 0;
    foreach my $file (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      $file =~ s{//}{/}g;
      $self->{'created_files'}->{$file} ||= $save_from_deletion;
    }
  }

  my @files = keys %{$self->{'created_files'}};
  return \@files;
}

sub get_save_status {
  my ($self, $file) = @_;
  return exists $self->{'created_files'}->{$file} ? $self->{'created_files'}->{$file} : 0;
}

=head2 delete_files

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Function  : remove files and directories stored in files_to_delete
  Returntype: n/a
  Exceptions: throws if deletion was not successful
  Example   : $self->delete_files;

=cut

sub delete_files {
    my $self = shift;

    return if ($self->save_files_for_deletion);

    my $input_files = $self->input_files;
    my $output_files = $self->output_files;

    FILE:
    foreach my $file (@{$self->created_files}) {
      next FILE if (!$self->_running && $self->get_save_status($file));
      next FILE if (grep {$_ eq $file} @$input_files);
      next FILE if (!$self->_running && grep {$_ eq $file} @$output_files);
      if ( -d $file ) {
          delete_directory($file, 1);
      }
      elsif ( -f $file ) {
          delete_file($file, 1);
      }
    }
    return;
}

sub _running {
  my $self = shift;
  if (@_) {
    $self->{'_running'} = (shift) ? 1 : 0;
  }
  return $self->{'_running'};
}

sub _temp_dir {
  my $self = shift;
  if (@_) {
    $self->{'_temp_dir'} = shift;
  }
  return $self->{'_temp_dir'};
}

sub get_temp_dir {
    my $self = shift;
    my $temp_dir = $self->_temp_dir;
    if (! $temp_dir) {
      $temp_dir = $self->working_dir()
            .'/'.$self->job_name.'.'.$$.'.tmp/';
      check_file_does_not_exist($temp_dir);
      $self->created_files($temp_dir);
      $self->_temp_dir($temp_dir);
      make_directory($temp_dir);
    }
    return $temp_dir;
}


sub options {
    my ($self, $option_name, $option_value) = @_;

    throw( "option_name not specified")
        if (! $option_name);

    $self->{'options'} ||= {};
    if (defined $option_value) {
        $self->{'options'}->{$option_name} = $option_value;
    }

    return $self->{'options'}->{$option_name};
}



1;
