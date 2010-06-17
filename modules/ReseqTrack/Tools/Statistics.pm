=pod

=head1 NAME

ReseqTrack::Tools::RunAlignment

=head1 SYNOPSIS

This is a base class for RunAlignment objects and provides some standard
accessor methods and throws exceptions when vital methods aren't implemented in
the child classes. The Child classes should wrap specific alignment algorithms

=head1 Example


=cut

package ReseqTrack::Tools::Statistics;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils;

sub new {
  my ($class, @args) = @_;
  my $self ={};
  bless $self,$class;
  my ($db, $skip_p2, $skip_p3) = 
    rearrange([qw(DB NEW_INDEX OLD_INDEX SKIP_P2 SKIP_P3)], @args);

  # setting defaults
  $skip_p2 = 1 unless(defined($skip_p2));
  $skip_p3 = 1 unless(defined($skip_p3));
  #####

  $self->db($db);
  $self->skip_p2($skip_p2);
  $self->skip_p3($skip_p3);

  return $self;
}

=head2 db

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Arg [2]   : ReseqTrack::DBSQL::DBAdaptor
  Function  : accessor method for the dbadaptor
  Returntype: ReseqTrack::DBSQL::DBAdaptor
  Exceptions: throws if not given a  ReseqTrack::DBSQL::DBAdaptor object
  Example   :

=cut



sub db{
  my ($self, $db) = @_;
  if($db){
    throw("Must pass SequenceIndexStatistics:db a ReseqTrack::DBSQL::DBAdaptor not ".
	  " $db") unless($db->isa("ReseqTrack::DBSQL::DBAdaptor"));
    $self->{db} = $db;
  }
  return $self->{db};
}



=head2 skip_p2/p3

  Arg [1]   : ReseqTrack::Tools::Statistics::SequenceIndexStatistics
  Arg [2]   : 0/1
  Function  : accessor method for skip_p2/p3, These flags indicate that pilot2
  or pilot3 runs should be skipped
  Returntype: 0/1
  Exceptions: n/a
  Example   : 

=cut



sub skip_p2{
  my ($self, $skip_p2) = @_;
  if($skip_p2){
    $self->{skip_p2} = $skip_p2;
  }
  return $self->{skip_p2};
}

sub skip_p3{
  my ($self, $skip_p3) = @_;
  if($skip_p3){
    $self->{skip_p3} = $skip_p3;
  }
  return $self->{skip_p3};
}

1;