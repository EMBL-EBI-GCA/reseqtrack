=pod

=head1 NAME

ReseqTrack::Tools::BatchSubmission

=head1 SYNOPSIS

This is a base class for BatchSubmission objects and provides some standard
accessor methods and throws exceptions when vital methods aren't implemented in
the child classes

=head1 Example


=cut

package ReseqTrack::Tools::BatchSubmission;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);


=head2 new

  Arg [1]   : ReseqTrack::Tools::BatchSubmission
  Arg [2]   : -program, string name of submission program
  Arg [3]   : -options, string commandline options for submission program
  Arg [4]   : -cmd, string, cmd submssion program should run
  Arg [5]   : -max_job_number, the maximum number of jobs this user can have in the
              system
  Arg [6]   : -user, string the username
  Function  : create a BatchSubmission object
  Returntype: ReseqTrack::Tools::BatchSubmission
  Exceptions: n/a
  Example   : 

=cut



sub new {
  my ($class, @args) = @_;
  my $self ={};
  bless $self,$class;
  my ($program, $options, $cmd, $max_job_number,
      $sleep, $user) = rearrange([qw(PROGRAM OPTIONS CMD MAX_JOB_NUMBER SLEEP USER)],
                                 @args);
  $self->program($program);
  $self->options($options);
  $self->cmd($cmd);
  $self->max_job_number($max_job_number);
  $self->sleep($sleep);
  $self->user($user);
  return $self;
}



=head2 program

  Arg [1]   : ReseqTrack::Tools::BatchSubmission
  Arg [2]   : string, name/path to executable submission program
  Function  : accessor method for program string
  Returntype: string
  Exceptions:
  Example   : 

=cut


sub program{
  my ($self, $program) = @_;
  if($program){
    $self->{'program'} = $program;
  }
  return $self->{'program'};
}



=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::BatchSubmission
  Arg [2]   : string, value for accessor
  Function  : These are accessor methods for the string based variables which are
  part of BatchSubmission, option, cmd, max_job_number, sleep and user. Options will 
  be any commandline options for the submission systems, cmd should be the command 
  to be run, max_job_number is the maximum number of jobs in the system at one time
  for this user, user defined the user and sleep is the amount of time the program
  sleeps for when the maximum number of jobs is reached
  Returntype: string
  Exceptions: 
  Example   : 

=cut



sub options{
  my ($self, $options) = @_;
  if($options){
    $self->{'options'} = $options;
  }
  return $self->{'options'};
}

sub cmd{
  my ($self, $cmd) = @_;
  if($cmd){
    $self->{'cmd'} = $cmd;
  }
  return $self->{'cmd'};
}

sub max_job_number{
  my ($self, $max_job_number) = @_;
  if($max_job_number){
    $self->{'max_job_number'} = $max_job_number;
  }
  return $self->{'max_job_number'};
}

sub sleep{
  my ($self, $sleep) = @_;
  if($sleep){
    $self->{'sleep'} = $sleep;
  }
  return $self->{'sleep'};
}

sub user{
  my ($self, $user) = @_;
  if($user){
    $self->{'user'} = $user;
  }
  return $self->{'user'};
}

=head2 construct_command_line

  Arg [1]   : ReseqTrack::Tools::BatchSubmission
  Function  : this method should construct a batch submission commandline.
  It should be implemented by the child class
  Returntype: string
  Exceptions: 
  Example   : my $cmd = construct_command_line();

=cut


sub construct_command_line{
  my ($self) = @_;
  throw($self." must implement a construct command line method");
}

=head2 get_total_job_number

  Arg [1]   : ReseqTrack::Tools::BatchSubmission
  Function  : This should calcualte the number of jobs this user has pending/running 
    on the system. It should be implemented by the child class
  Returntype: string
  Exceptions: 
  Example   : my $cmd = get_total_job_number();

=cut

sub get_total_job_number{
  my ($self) = @_;
  throw($self." must implement get total job number");
}

=head2 check_for_awol_jobs

  Arg [1]   : ReseqTrack::Tools::BatchSubmission
  Function  : This should check if any jobs have disappeared from the submission
    system without being marked as failed or successful. It should be implemented 
    by the child class
  Returntype: string
  Exceptions: 
  Example   : my $cmd = check_for_awol_jobs();

=cut

sub check_for_awol_jobs{
  my ($self) = @_;
  throw($self." must implement check for awol jobs");
}


=head2 has_to_many_jobs

  Arg [1]   : ReseqTrack::Tools::BatchSubmission
  Arg [2]   : int, maximum number of jobs
  Arg [3]   : int, seconds to sleep for
  Arg [4]   : string, username
  Function  : calculate how many jobs are running and only continue if the number
  is less than maximum after sleeping for an appropriate amount of time
  Returntype: n/a
  Exceptions: n/a
  Example   : 

=cut



sub has_to_many_jobs{
  my ($self, $max, $sleep, $user) = @_;
  $user = $self->user unless($user);
  my $total = get_total_job_number($user);
  $max = $self->max_job_number unless($max);
  $sleep = $self->sleep unless($sleep);
  if($total < $max){
    return 0;
  }else{
    print STDERR "Have ".$total." jobs compared to ".$max." going to sleep for ".
        $sleep." seconds\n";
    CORE::sleep($sleep);
    return has_to_many_jobs($max, $sleep, $user);
  }
}

1;
