package ReseqTrack::Tools::Loader::LoadFiles;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Loader);

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  return $self;
}
