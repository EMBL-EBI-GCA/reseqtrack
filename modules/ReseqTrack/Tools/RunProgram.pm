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


my %files_to_delete;
my @command_history;
my $num_run_programs = 0;

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

sub new {
  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;

  #input files, we will allow single file name or list of file names
  my ( $input_files,  $program, $job_name, $working_dir,
        $echo_cmd_line, $save_files_for_deletion )
    = rearrange( [
         qw( INPUT_FILES PROGRAM JOB_NAME WORKING_DIR
             ECHO_CMD_LINE SAVE_FILES_FOR_DELETION )
		], @args);

  if ( ! $working_dir){
    $working_dir = "/tmp/";
  }

  $self->input_files($input_files);
  $self->program($program);
  $self->job_name($job_name);
  $self->working_dir($working_dir);
  $self->echo_cmd_line($echo_cmd_line);
  $self->save_files_for_deletion($save_files_for_deletion);

  $num_run_programs++;

  return $self;
}

=head2 DESTROY

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Function  : calls delete_files() if it is the last RunProgram object to exist
  Returntype: 
  Exceptions: 
  Example   : 

=cut

sub DESTROY {
    my $self = shift;

    $num_run_programs--;

    if ( $num_run_programs == 0) {
        if (! $self->save_files_for_deletion ) {
            $self->delete_files();
        }
    }
}

=head2 run

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Function  : each child object should implement a run method
  Returntype: 
  Exceptions: throws as this method should be implemented in the child class
  Example   : 

=cut

sub run {
  my ($self) = @_;
  throw(  $self
          . " must implement a run method as ReseqTrack::Tools::RunProgram "
          . "does not provide one" );
}

=head2 execute_command_line

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, command line to execute
  Function  : executes command, stores command in command_history
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

    $self->command_history($command_line);

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
    $self->{'save_files_for_deletion'} = (shift) ? 1 : 0;
  }
  return $self->{'save_files_for_deletion'};
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
    $self->{'echo_cmd_line'} = (shift) ? 1 : 0;
  }
  return $self->{'echo_cmd_line'};
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
    check_file_exists($arg);
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

=head2 command_history

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, executed command
  Function  : accessor method for the history of commands executed by all RunProgram objects
  Returntype: arrayref of strings
  Exceptions: 
  Example   : my $first_command = ${$self->command_history}[0];

=cut

sub command_history {
    my ($self, $command_line) = @_;

    if (defined $command_line) {
        push(@command_history, $command_line);
    }
    return \@command_history;
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

    my $job_name = $$;

    $self->job_name($job_name);
    return $job_name;
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
    $self->{'working_dir'} = $arg;
  }
  return $self->{'working_dir'};
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

  if (! -d $dir) {
      make_directory($dir);
  }

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

  if (! $self->{'output_files'}) {
      $self->{'output_files'} = [] ;
  }

  if ($arg) {
    if ( ref($arg) eq 'ARRAY' ) {
      push( @{ $self->{'output_files'} }, @$arg );
    } else {
      push( @{ $self->{'output_files'} }, $arg );
    }
  }
  return $self->{'output_files'};
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

  if (! $self->{'input_files'}) {
      $self->{'input_files'} = [] ;
  }

  if ($arg) {
    if ( ref($arg) eq 'ARRAY' ) {
        foreach my $file (@$arg) {
            check_file_exists($file);
            push( @{ $self->{'input_files'} }, $file );
        }
    } else {
        check_file_exists($arg);
        push( @{ $self->{'input_files'} }, $arg );
    }
  }

  return $self->{'input_files'};
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

sub files_to_delete {
  my ( $self, $file ) = @_;

  if ($file) {
    if ( ref($file) eq 'ARRAY' ) {
      foreach my $path (@$file) {
	$files_to_delete{$path} = 1;
      }
    } else {
      $files_to_delete{$file} = 1;
    }
  }
  my @keys = keys( %files_to_delete );
  return \@keys;
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

    foreach my $file (@{$self->files_to_delete}) {
        print "Deleting " . $file . "\n";

        if ( -d $file ) {
            delete_directory($file);
        }
        else {
            delete_file($file);
        }

        delete $files_to_delete{$file};
    }

    return;
}

1;
