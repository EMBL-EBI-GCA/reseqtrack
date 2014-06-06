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

use ReseqTrack::Tools::FileSystemUtils
    qw(check_file_exists check_directory_exists delete_directory
    delete_file check_file_does_not_exist check_executable create_tmp_process_dir);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);
use File::Find qw(find);
use Env qw( @PATH );
use POSIX;


=head2 new

  Arg [-input_files]   :
      string or arrayref of strings,  paths to any input files
  Arg [-program]   :
      string, the program (if in $PATH) or path to the program
  Arg [-job_name]   :
      string, name of the job. This will be generated automatically if not specified.
  Arg [-working_dir]   :
      string, default '/tmp/', path of the working directory
  Arg [-options]   :
      hashref, {option_name => option_value}
      options specific to the child class
  Arg [-echo_cmd_line]   :
      boolean flag, default 0, echo the command line rather than execute commands
      this option is designed for debugging purposes
  Arg [-save_files_from_deletion]   :
      boolean flag, default 0, save temporary files (i.e. do not delete them)
      this option is designed for debugging purposes

  Arg [-output_name_prefix] :
  	  string, default is "RESULTS"; it is prefix to be put in output files
  	  
  Function  : Creates a new ReseqTrack::Tools::RunProgram object 
  Returntype: ReseqTrack::Tools::RunProgram
  Exceptions: 
  Example   : my $run_program = ReseqTrack::Tools::RunProgram->new(
                -input_files => ['/path/file1', '/path/file2'],
                -program => "program",
                -working_dir => '/path/to/dir/',
                -job_name => "my_job",
                -options => {"flag1" => 1, "flag2" => 0}
                -output_name_prefix => 'PHASE1'
                );

=cut


# package variables:
my $GLOBAL_ECHO_CMD_LINE;
my $GLOBAL_save_files_from_deletion;
my $GLOBAL_WORKING_DIR;

my $term_sig =  0;
$SIG{TERM} = \&termhandler;
$SIG{INT} = \&termhandler;
sub termhandler {
    $term_sig = 1;
}

sub new {
  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;

  #input files, we will allow single file name or list of file names
  my ( $input_files,  $program, $job_name, $working_dir,
        $options,
        $echo_cmd_line, $save_files_from_deletion, $output_name_prefix )
    = rearrange( [
         qw( INPUT_FILES PROGRAM JOB_NAME WORKING_DIR
             OPTIONS
             ECHO_CMD_LINE save_files_from_deletion OUTPUT_NAME_PREFIX)
		], @args);

  $self->input_files($input_files);
  $self->program($program);
  $self->job_name($job_name);
  $self->working_dir($working_dir);
  $self->echo_cmd_line($echo_cmd_line);
  $self->save_files_from_deletion($save_files_from_deletion);
  $self->output_name_prefix($output_name_prefix);

  $self->options($self->DEFAULT_OPTIONS);
  $self->options($options) if ($options);
  
  return $self;
}

=head2 DEFAULT_OPTIONS

  Function  : child classess can implement a DEFAULT_OPTIONS method which defines the
              default values for $self->options
  Returntype: hashref

=cut

sub DEFAULT_OPTIONS {return {};}

=head2 run

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : Any additional args will be passed to the child class
  Function  : Calls the run_program method, which must be implemented by the child class
              Does various checks that are needed for running any program
  Returntype: 
  Exceptions: Throws if program executable does not exist.
              Throws if any of the input files do not exist.
  Example   : $my_program_object->run;

=cut

sub run {
  my ($self, @args) = @_;

  $self->_running(1);

  throw("new policy for RunProgram: output directory must always be set explicity") if !$self->working_dir;

  my $program = $self->program;
  throw "do not have a program executable\n" if (! $program);
  if (! -d $program) {
    check_executable($program);
  }

  foreach my $file (@{$self->input_files}) {
    check_file_exists($file);
  }
  if (! $self->job_name) {
      $self->generate_job_name;
  }

  $self->change_dir;
  my @returned_values = $self->run_program(@args);

  $self->_running(0);
  $self->delete_files;

  return @returned_values;
}

=head2 DESTROY

  Function  : ensures that files will be deleted when a RunProgram object goes out of scope
              This can happen either after successful completion or if the perl script quits 
              unexpectedly e.g. when an exception is thrown.

=cut

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
    my ($self, $command_line, $exit_code_handler) = @_;

    $exit_code_handler //= \&default_process_exit_code;

    if ($self->echo_cmd_line) {
        $command_line = "echo \'" . $command_line . "\'";
    }
    print "Executing command:$/";
    print $command_line . $/;

    my $pid = fork;
    if (!defined $pid) {
      throw "could not fork: $!";
    }
    elsif(!$pid) {
      setpgrp(0,0);
      $self->save_files_from_deletion(1);
      exec(bash => (-o => 'pipefail', -c => $command_line)) or do{
        print STDERR "Did not execute: $!: $command_line\n";
        exit(255);
      }
    }
    while (! waitpid($pid, WNOHANG)) {
      sleep(3);
      if ($term_sig) {
        $self->term_sig($term_sig);
        my $deleting_pid = fork;
        if ($deleting_pid == 0) {
          kill -15, $pid;
          #`kill -s 15 -$pid`;
          $self->delete_files;
          exit(0);
        }
        $self->save_files_from_deletion(1);
        throw("received a signal ($term_sig) so command line was killed $command_line");
      }
    }
    #throw("received a signal ($term_sig) so command line was killed $command_line") if ($term_sig);

    my $signal = $? & 127;
    throw("process died with signal $signal $command_line") if ($signal);
    my $exit = $? >> 8;
    throw("command could not be executed by bash: $command_line") if ($exit == 255);
    $exit_code_handler->($self, $exit, $command_line);
    return $exit;
}

sub default_process_exit_code {
  my ($self, $exit, $command_line) = @_;
  throw("command exited with value $exit $command_line") if ($exit != 0);
}

=head2 save_files_from_deletion

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : boolean, optional, value of save_files_from_deletion flag
  Function  : accessor method for save_files_from_deletion flag.  Files marked for deletion
              will not be deleted. Designed to be used for debugging purposes.
  Returntype: boolean
  Exceptions: 
  Example   : my $flag = $self->save_files_from_deletion;

=cut

sub save_files_from_deletion {
  my $self = shift;
  if (@_) {
    $GLOBAL_save_files_from_deletion = (shift) ? 1 : 0;
  }
  return $GLOBAL_save_files_from_deletion;
}

=head2 echo_cmd_line

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : boolean, optional, value of echo_cmd_line flag
  Function  : accessor method for echo_cmd_line flag.  Command lines will be 'echoed'
              but not executed. Designed to be used for debugging purposes.
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


=head2 output_dir

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, path of output directory
  Function  : alternative accessor to working_dir.  output_dir and working_dir are the same thing.
  Returntype: string
  Exceptions: n/a
  Example   : my $output_dir = $self->output_dir;

=cut

sub output_dir {
  my ( $self, $arg ) = @_;
  return $self->working_dir($arg);
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


=head2 output_name_prefix

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, optional, a prefix to be put in an output file name to identify the results of a particular run
  Function  : accessor method for output_name_prefix
  Returntype: strins
  Exceptions: n/a
  Example   : $self->output_name_prefix('PHASE1');

=cut

sub output_name_prefix {
  my ($self, $output_name_prefix) = @_;
  if ($output_name_prefix) {
    $self->{'output_name_prefix'} = $output_name_prefix;
  }
  return $self->{'output_name_prefix'};
}

=head2 output_files

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for output files
              output files are also added to the list of created files but marked to be saved from deletion
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
  Function  : accessor method for input files.
  Returntype: arrayref of strings
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

=head2 created_files

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string or arrayref of strings
  Arg [3]   : boolean, save_file_from_deletion, default 0
  Function  : Store files and directories that have been created by the program or class
              These files will be deleted when the object is destroyed.
              Files marked with the save_file_from_deletion flag will only be saved if the run method
              completes successfully i.e. they will be deleted if an error occurs
  Returntype: arrayref of strings
  Example   : $self->created_files('path/to/file');

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

=head2 get_save_status

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, filename
  Function  : returns true if a file has been marked to be saved from deletion.
              This is called by the delete_files method.
  Returntype: boolean
  Example   : $save_file = $self->get_save_status('path/to/file');

=cut

sub get_save_status {
  my ($self, $file) = @_;
  return exists $self->{'created_files'}->{$file} ? $self->{'created_files'}->{$file} : 0;
}

=head2 delete_files

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Function  : remove files and directories stored in created_files
              Files marked as save_from_deletion are also deleted if the program quits unexpectedly.
  Returntype: n/a
  Exceptions: throws if deletion was not successful
  Example   : $self->delete_files;

=cut

sub delete_files {
    my $self = shift;

    return if ($self->save_files_from_deletion);

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

=head2 _running

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : boolean, program is running
  Function  : accessor method for the _running flag.
              This is set and unset by the run method
              ALL files will be deleted if the RunProgram object is destroyed while this flag is set to 1.
  Returntype: boolean
  Example   : my $flag = $self->_running;

=cut

sub _running {
  my $self = shift;
  if (@_) {
    $self->{'_running'} = (shift) ? 1 : 0;
  }
  return $self->{'_running'};
}

=head2 get_temp_dir

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Function  : Gets a temp directory for use by the child class.
              Adds the temp directory to list of created files.
  Returntype: String, path of temp directory
  Exceptions: Throws if there is an error making the temp directory.
  Example   : my $dir = $self->get_temp_dir;

=cut

sub get_temp_dir {
    my ($self) = @_;
    my $temp_dir = $self->{'_temp_dir'};
    if (! $temp_dir) {
      $temp_dir = create_tmp_process_dir($self->working_dir, $self->job_name, 0);
      $self->created_files($temp_dir);
      $self->{'_temp_dir'} = $temp_dir;
    }
    return $temp_dir;
}

=head2 get_short_input_names

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [1]   : 1 or 2 to indicate to which level of shortness the name should be; with 2, the file name will be just an index.  Default is 1. 
  Function  : Uses symbolic links to allow the conversion of long file names into something shorter
              Symbolic links are created within a temporary directory
  Returntype: hash, key is the long file name, value is the short file name
  Exceptions: Throws if there is an error linking the files
  Example   : my $short_name = $self->get_short_input_names->{'/long/file/name'};
	      OR
              my $short_names = $self->get_short_input_names(2);
	      my $short_name = $short_names->{'/long/file/name'};

=cut


sub get_short_input_names {
    my ($self, $level) = @_;
    $level = 1 if (!$level);
    my $temp_dir = $self->get_temp_dir;
    my $temp_dir_basename = fileparse($temp_dir);
    $self->{'_short_inputs'} ||= {};
    my $cnt = 0;    
    foreach my $long_name (grep {!$self->{'_short_inputs'}->{$_}} @{$self->input_files}) {
      my ($basename, $dirname) = fileparse($long_name);
      $dirname =~ s/\/$//;
      my $short_name;
      if ( $level == 1) { 
      	$short_name = $temp_dir_basename . '/' . $basename;
      }
      elsif ($level == 2) {
	 $short_name = $temp_dir_basename . '/' . $cnt;
      }
      else {
	throw("Please use either 1 or 2 as level");
      }
      $self->{'_short_inputs'}->{$long_name} = $short_name;
      symlink($long_name, $short_name) or throw("could not symlink $long_name to $short_name");

      my $link_others_sub = sub {
        return if $File::Find::dir ne $dirname;
        my ($extension) = $_ =~ /$basename(.+)$/;
        return if ! defined $extension;
        my $destination = $temp_dir . '/' . $_;
        #my $destination = $temp_dir . '/' . $cnt;
	symlink($File::Find::name, $destination) or throw("could not symlink $File::Find::name to $destination $!");
      };
      find($link_others_sub, $dirname);
      $cnt++;
    }
    return $self->{'_short_inputs'};
}


=head2 options

  Arg [1]   : ReseqTrack::Tools::RunProgram
  Arg [2]   : string, option_name, or a hash ref of options
  Arg [3]   : any, option_value (optional)
  Function  : Accessor method for options required by the child class
  Returntype: option_value, or hash ref of all options
  Exceptions: Throws if option_name is not specified.
  Example   : my $option_value = $self->options('option_name');
  Example   : my $option_value = $self->options('option_name','option_value');
  Example	: my $options_hashref = $self->options();
  Example   : my $options_hashref = $self->options(\%options); # will merge these options with any already set
  
=cut


sub options {
  my ($self, @args) = @_;

  $self->{'options'} ||= {};  
  my $num_of_args = scalar(@args);

  # no arguments, return all the options
  if ($num_of_args == 0){
    return $self->{'options'};
  }
  elsif ($num_of_args == 1){
    my $ref_type = ref($args[0]);

    if (! $ref_type){
      # arg is a scalar
      return $self->{'options'}->{$args[0]};
    }
    elsif ($ref_type eq 'HASH'){
      # merge these options with any existing

      while (my ($name, $value) = each %{$args[0]}) {
        $self->{'options'}->{$name} = $value;
      }

      return $self->{'options'};
    }
    else {
      throw("Cannot set options with a $ref_type reference");
    }

  }
  else{
    my ($option_name,$option_value) = @args;
    throw( "option_name not specified") if (! $option_name);
    $self->{'options'}->{$option_name} = $option_value;
    return $self->{'options'}->{$option_name};
  }

}

sub term_sig {
  my ($self, $term_sig) = @_;
  if (@_ > 1) {
    $self->{'_term_sig'} = $term_sig
  }
  return $self->{'_term_sig'} || 0;
}


1;
