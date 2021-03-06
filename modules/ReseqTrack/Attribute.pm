package ReseqTrack::Attribute;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use Data::Dumper;
use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($table_name, $other_id, $attribute_name, $attribute_value, $attribute_units ) =
      rearrange([qw(TABLE_NAME OTHER_ID ATTRIBUTE_NAME ATTRIBUTE_VALUE ATTRIBUTE_UNITS) ], @args);
  $self->table_name($table_name);
  $self->other_id($other_id);
  $self->attribute_name($attribute_name);
  $self->attribute_value($attribute_value);
  $self->attribute_units($attribute_units);
  return $self;
}



=head2 accessor methods

  Arg [1]   : Reseq::Statistics
  Arg [2]   : mostly strings
  Function  : set variable in object
  Returntype: return variable
  Exceptions: n/a
  Example   : 

=cut

sub other_id{
  my ($self, $arg) = @_;
  if($arg){
    $self->{other_id} = $arg;
  }
  return $self->{other_id};
}

sub table_name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{table_name} = $arg;
  }
  return $self->{table_name};
}

sub attribute_name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{attribute_name} = $arg;
  }
  return $self->{attribute_name};
}

sub attribute_value{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{attribute_value} = $arg;
  }
  return $self->{attribute_value};
}

sub attribute_units{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{attribute_units} = $arg;
  }
  return $self->{attribute_units};
}





1;
