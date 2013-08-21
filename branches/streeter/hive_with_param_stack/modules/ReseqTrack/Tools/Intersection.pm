
=head1 NAME

ReseqTrack::Tools::Intersection

=head1 SYNOPSIS

     my $list_set = ReseqTrack::Tools::Intersection
                     ->new(
                           -LIST => \@list,
                          );
     my $unique = $list_set->not($second_set);

=head1 DESCRIPTION

This is a object which will hold a list of strings and compare them
using and, not, or and xor

=head1 CONTACT

 laura@ebi.ac.uk or zi@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details exported static class methods. 

=cut


package ReseqTrack::Tools::Intersection;

use vars qw(@ISA);
use strict;
use Carp;
use warnings;
use ReseqTrack::Tools::Argument qw( rearrange );

@ISA = qw();


=head2 new

  Arg [1]   : LIST, listref
  Function  : instantiate a new ReseqTrack::Tools::Intersection object
  Returntype: ReseqTrack::Tools::Intersection
  Exceptions: if the value its passed for LIST isn't a listref'
  Caller    :
  Example   : ReseqTrack::Tools::Intersection->new();

=cut


sub new{
  my ($class, @args) = @_;
  my $self = bless {}, $class;
 
  #print "Args = ".join("\t", @args)."\n";
  my ($list) = rearrange([qw(LIST)], @args);

  $self->{'_ID_list'} = [];

  if($list){
    if(ref($list) eq 'ARRAY'){
      $self->{'_ID_list'} = $list;
    }else{
      confess("LIST need to be a listref not a $list");
    }
  }
  return $self;
}


=head2 add/delete_element

  Arg [1]   : array
  Function  : adds to the list or removes from the list a set of elements
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : $self->add_element('contig_id');

=cut



sub add_element{
  my ($self, @elements) = @_;
  if(@elements){
    push(@{$self->{'_ID_list'}}, @elements);
  }
}


sub delete_element{
  my ($self, @elements) = @_;
  my @ids = @{$self->{'_ID_list'}};
  foreach my $e(@elements){
    my @new_ids = grep $_ != $e, @ids;
    @ids = @new_ids;
    @new_ids = ();
  }

  @{$self->{'_ID_list'}} = @ids;
}


=head2 list

  Arg [1]   : none
  Function  :
  Returntype: returns the listref to the caller
  Exceptions: none
  Caller    :
  Example   : my @array = @{$idset->ID_list};

=cut



sub list{
  my ($self) = @_;
  return $self->{'_ID_list'};
}



=head2 _calculate

  Arg [1]   : ReseqTrack::Tools::Intersection
  Function  : to produce a hash which can be used for unions and intersections
  Returntype: ReseqTrack::Tools::Intersection
  Exceptions: if the value passed in isn't a ReseqTrack::Tools::Intersection'
  Caller    : 
  Example   : $self->_calculate($idlist);

=cut



sub _calculate{
  my ($self, $idlist) = @_;

  if($idlist){
    $idlist->isa('ReseqTrack::Tools::Intersection') ||
    confess("Need to pass these methods (and/or/not) an " .
		 "ReseqTrack::Tools::Intersection object otherwise won't work " .
		 "you passed in $idlist");
  }else{
    confess("you must pass in an ReseqTrack::Tools::Intersection object to ".
          "the and/or method not".
          "$idlist an undefined value");
  }

  my %count;
  my @own_ids =  @{$self->{'_ID_list'}};
  my @comp_ids = @{$idlist->list};
  foreach my $e(@own_ids, @comp_ids) {$count{$e}++}
  return \%count;
}


=head2 and

  Arg [1]   : ReseqTrack::Tools::Intersection
  Function  : to produce a list which contains every id which exists 
  on both lists. This is done by using calculate to produce a hash which 
  contains the elements on both lists together with the number of times each
  element appears if an elements appears twice it must live on both lists so
  is added to the new list (set intersection)
  Returntype: ReseqTrack::Tools::Intersection
  Exceptions: none
  Caller    : 
  Example   : my $c = $a->and($b);

=cut



sub and{
  my ($self, $idlist) = @_;
 
  my %count = %{$self->_calculate($idlist)};
  my @and;
  foreach my $e(keys(%count)){
    if($count{$e} == 2){
      push(@and, $e);
    }
  }
  my $idset = ReseqTrack::Tools::Intersection->new(
						 -LIST => \@and,
						);
  return $idset;
}


=head2 or

  Arg [1]   : ReseqTrack::Tools::Intersection
  Function  : to return all the values which appear in at least on of the 
  lists this is done by caling calculate and returning an idlist of all the
  keys of the hash it returns as this all live on at least one of the lists
  (set union)
  Returntype: ReseqTrack::Tools::Intersection
  Exceptions: none
  Caller    : 
  Example   : $c = $a->or($b);

=cut



sub or{
  my ($self, $idlist) = @_;
  
  my %count = %{$self->_calculate($idlist)};
  my @and = keys(%count);
  my ($p, $f, $l) = caller;
  my $idset = ReseqTrack::Tools::Intersection->new(
						 -LIST => \@and,
						);
  return $idset;
}


=head2 not

  Arg [1]   : ReseqTrack::Tools::Intersection
  Function  : to return all the ids on the objects current list but not on
  the list passed to the function this is the simple difference this is done
  by creating a hash of the list for comparision and returning an ReseqTrack::Tools::Intersection of
  all the elements in the current IDset which aren't in that hash'
  Returntype: ReseqTrack::Tools::Intersection
  Exceptions: confesss if not passed an ReseqTrack::Tools::Intersection
  Caller    :
  Example   : $c = $a->not($b);

=cut


#as a not given all these methods return ReseqTrack::Tools::Intersections they can be chained together
#to return a set of things present in two list or not a third for example
# $d = $a->not($b)->and($c);

sub not{
  my ($self, $idlist) = @_;
  if($idlist){
    $idlist->isa('ReseqTrack::Tools::Intersection') ||
    confess("Need to passs these methods (and/or/not) an " .
		 "ReseqTrack::Tools::Intersection object otherwise won't work " .
		 "you passed in $idlist");
  }else{
    confess("you must pass in an ReseqTrack::Tools::Intersection object to the not method not".
		 "$idlist an undefined value");
  }
  my @own_ids =  @{$self->{'_ID_list'}};
  my @comp_ids = @{$idlist->list};
  my %seen;
  my @and;
  foreach my $e(@comp_ids) {
    $seen{$e}++;
  }
  foreach my $e(@own_ids){
    if(exists($seen{$e})){
      next;
    }
    push(@and, $e)
  }
  my $idset = ReseqTrack::Tools::Intersection->new(
						 -LIST => \@and,
						);
  return $idset;
}


=head2 xor

  Arg [1]   : ReseqTrack::Tools::Intersection
  Function  : produces a IDset which contains all the IDs which exist on only
  one list or the other and not both (symmetric difference)
  Returntype: ReseqTrack::Tools::Intersection
  Exceptions: none
  Caller    :
  Example   : $c = $a->xor($b);

=cut



sub xor {
  my ($self, $idlist) = @_;

  my %count = %{$self->_calculate($idlist)};
  my @and;
  foreach my $e(keys(%count)){
    if($count{$e} == 1){
      push(@and, $e);
    }
  }
  my $idset = ReseqTrack::Tools::Intersection->new(
						 -LIST => \@and,
						);
  return $idset;
}


=head2 count

  Arg [1]   : none
  Function  : returns length of list
  Returntype: integer
  Exceptions: none
  Caller    : 
  Example   : my $count = $idset->count;

=cut


sub count{
  my ($self) = @_;

  if(@{$self->{'_ID_list'}}){
     return scalar(@{$self->{'_ID_list'}});
  }else{
    return 0;
  }
}


=head2 subset

  Arg [1]   : size, int
  Function  : produces a slice/subset of the list
  Returntype: ReseqTrack::Tools::Intersection
  Exceptions: confesss if not passes a size
  Caller    : 
  Example   : my $subset = $idset->subset(10);

=cut



sub subset{
  my ($self, $size) = @_;

  confess("need a size inorder to return a subset of $self") 
    unless $size;
  if($size > $self->count){
    $size = $self->count;
  };
  my @array = @{$self->list};
  my $end = $size -1;
  my @subset = @array[0..$end];

  my $idset = ReseqTrack::Tools::Intersection->new(
                         -LIST => \@subset,
                        );

  return $idset;
  
}



1;
