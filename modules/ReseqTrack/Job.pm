
=pod

=head1 NAME

ReseqTrack::Job

=head1 SYNOPSIS

A container to hold the data for the Job and Job status table which are used to
track jobs during in the event pipeline system

=head1 EXAMPLE

$job = ReseqTrack::Job->new
        (
         -input_string => $input_string,
         -event => $event,
         -retry_count => 0,
         -current_status => 'CREATED',
        );

=cut

package ReseqTrack::Job;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Base;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::Base);

=head2 new

  Arg [1]   : ReseqTrack::Job
  Arg [2]   : int, submission id from the batch submission system used
  Arg [3]   : ReseqTrack::Event object
  Arg [4]   : string, input string
  Arg [5]   : string, stdout filepath
  Arg [6]   : string, stderr filepath
  Arg [8]   : string, the machine the job is run on
  Arg [9]   : int, retry count, how many times the job has been attempted
  Arg [10]  : string, current status, the current status of the job, this is stored
  in the job status table and can be CREATED, SUBMITTED, WAITING, RUNNING, FAILED
  AWOL or SUCCESSFUL
  Function  : create ReseqTrack::Job object
  Returntype: ReseqTrack::Job
  Exceptions: throws if not given an event object or an input string
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my ( $p, $f, $l ) = caller;

    my (
        $submission_id, $event,          $input,
        $stdout,        $stderr,         $host,
        $retry_count,   $current_status, $time
      )
      = rearrange(
        [
            qw( SUBMISSION_ID EVENT INPUT_STRING
              STDOUT STDERR EXEC_HOST RETRY_COUNT CURRENT_STATUS TIME)
        ],
        @args
      );

    #Error handling
    throw("Can't create a job without an event object") unless ($event);
    throw("Can't create a job without a input string")  unless ($input);
    #####

    $self->{input_string}   = $input;
    $self->{submission_id}  = $submission_id;
    $self->{event}          = $event;
    $self->{stdout_file}    = $stdout;
    $self->{stderr_file}    = $stderr;
    $self->{host}           = $host;
    $self->{retry_count}    = $retry_count;
    $self->{current_status} = $current_status;
    $self->{time}           = $time;
    return $self;
}

=head2 accessor methods

  Arg [1]   : Reseq::History
  Arg [2]   : mostly strings
  Function  : set variable in object
  Returntype: return variable
  Exceptions: n/a
  Example   : 

=cut

sub submission_id {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{submission_id} = $arg;
    }
    return $self->{submission_id};
}

sub event {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{event} = $arg;
    }
    return $self->{event};
}

sub input_string {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{input_string} = $arg;
    }
    return $self->{input_string};
}

sub time {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{time} = $arg;
    }
    return $self->{time};
}

sub host {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{host} = $arg;
    }
    return $self->{host};
}

sub retry_count {
    my ( $self, $arg ) = @_;
    if ( defined $arg ) {
        $self->{retry_count} = $arg;
    }
    return $self->{retry_count};
}

=head2 stdout/err_file

  Arg [1]   : ReseqTrack::Job
  Arg [2]   : string, path to stdout or err file
  Function  : set and return the string for the stdout or err files. It will also
  remove any // from the paths
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub stdout_file {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $arg =~ s/\/\//\//g;
        $self->{stdout_file} = $arg;
    }
    return $self->{stdout_file};
}

sub stderr_file {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $arg =~ s/\/\//\//g;
        $self->{stderr_file} = $arg;
    }
    return $self->{stderr_file};
}

=head2 current_status

  Arg [1]   : ReseqTrack::Job
  Arg [2]   : string, status
  Arg [3]   : binary flag, 0/1 to indicate if the status should be updated in the 
  database
  Function  : store the give status or retrieve the status from the database if
   not defined or the refetch flag is set
  Returntype: string
  Exceptions: n/a
  Example   : 

=cut

sub current_status {
    my ( $self, $arg, $refetch ) = @_;
    if ($arg) {
        $self->{current_status} = $arg;
    }
    if ($refetch) {
        $self->{current_status} = undef;
    }
    if ( $self->adaptor && !$self->{current_status} ) {
        $self->adaptor->fetch_status($self);
    }
    return $self->{current_status};
}

1;

