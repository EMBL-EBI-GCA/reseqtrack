package ReseqTrack::Job;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Base;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::Base);

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($p, $f, $l) = caller;
  #print "Creating job with ".join(" ", @args)." $f:$l\n";
  my ($submission_id, $event, $input, $stdout, $stderr, 
      $host, $retry_count, $current_status, $time) =
  rearrange([qw( SUBMISSION_ID EVENT INPUT_STRING 
                 STDOUT STDERR EXEC_HOST RETRY_COUNT CURRENT_STATUS TIME) ], @args);
  #Error handling
  throw("Can't create a job without an event object") unless($event);
  throw("Can't create a job without a input string") unless($input);
  #####

  $self->{input_string} = $input;
  $self->{submission_id} = $submission_id;
  $self->{event} = $event;
  $self->{stdout_file}  = $stdout;
  $self->{stderr_file} = $stderr;
  $self->{host} = $host;
  $self->{retry_count}  = $retry_count;
  $self->{current_status} = $current_status;
  $self->{time} = $time;
  return $self;
}

#submission id
#event
#array index
#input_string
#stdout file
#stderr file
#host
#retry count
#current_status
#time


sub submission_id{
  my ($self, $arg) = @_;
  if($arg){
    $self->{submission_id} = $arg;
  }
  return $self->{submission_id};
}
 
sub event{
  my ($self, $arg) = @_;
  if($arg){
    $self->{event} = $arg;
  }
  return $self->{event};
}


 sub input_string{
  my ($self, $arg) = @_;
  if($arg){
    $self->{input_string} = $arg;
  }
  return $self->{input_string};
}

 sub stdout_file{
  my ($self, $arg) = @_;
  if($arg){
    $arg =~ s/\/\//\//g ;
    $self->{stdout_file} = $arg;
  }
  return $self->{stdout_file};
}

 sub stderr_file{
  my ($self, $arg) = @_;
  if($arg){
    $arg =~ s/\/\//\//g ;
    $self->{stderr_file} = $arg;
  }
  return $self->{stderr_file};
}

 sub host{
  my ($self, $arg) = @_;
  if($arg){
    $self->{host} = $arg;
  }
  return $self->{host};
}
 
sub retry_count{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{retry_count} = $arg;
  }
  return $self->{retry_count};
}

 sub current_status{
  my ($self, $arg, $refetch) = @_;
  if($arg){
    $self->{current_status} = $arg;
  }
  if($refetch){
    $self->{current_status} = undef;
  }
  if($self->adaptor && !$self->{current_status}){
    $self->adaptor->fetch_status($self);
  }
  return $self->{current_status};
}

sub time{
  my ($self, $arg) = @_;
  if($arg){
    $self->{time} = $arg;
  }
  return $self->{time};
}
 
1;


