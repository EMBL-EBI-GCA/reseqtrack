package ReseqTrack::RejectLog;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($file_id, $is_reject, $reject_reason, $created ) = 
      rearrange([qw(FILE_ID IS_REJECT REJECT_REASON CREATED ) ], @args);
  $self->file_id($file_id);
  $self->is_reject($is_reject);
  $self->reject_reason($reject_reason);
  $self->created($created);
  return $self;
}

=head2 accessor methods

  Arg [1]   : Reseq::RejectLog
  Arg [2]   : mostly strings
  Function  : set variable in object
  Returntype: return variable
  Exceptions: n/a
  Example   : 

=cut

sub file_id{
  my ($self, $arg) = @_;
  if($arg){
    $self->{file_id} = $arg;
  }
  return $self->{file_id};
}

sub is_reject{
  my ($self, $arg) = @_;
  if($arg){
    $self->{is_reject} = $arg;
  }
  return $self->{is_reject};
}

sub reject_reason{
  my ($self, $arg) = @_;
  if($arg){
    $self->{reject_reason} = $arg;
  }
  return $self->{reject_reason};
}

sub created{
  my ($self, $arg) = @_;
  if($arg){
    $self->{created} = $arg;
  }
  return $self->{created};
}


1;
