=pod

=head1 NAME

ReseqTrack::Tools::Statistics

=head1 SYNOPSIS

Base class for generating statistics on the basis of sequence or alignment indexes

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
  my ($db, $new_index, $old_index, $skip_p2, $skip_p3, $collection_type, 
      $collection_name) = rearrange([qw(DB NEW_INDEX OLD_INDEX SKIP_P2 SKIP_P3 
					COLLECTION_TYPE	COLLECTION_NAME)], @args);

  # setting defaults
  $skip_p2 = 1 unless(defined($skip_p2));
  $skip_p3 = 1 unless(defined($skip_p3));
  #####

  $self->db($db);
  $self->new_index($new_index);
  $self->old_index($old_index);
  $self->skip_p2($skip_p2);
  $self->skip_p3($skip_p3);
  $self->collection_type($collection_type);
  $self->collection_name($collection_name);
  return $self;
}

=head2 new/old_index

  Arg [1]   : ReseqTrack::Tools::Statistics
  Arg [2]   : string, filepath
  Function  : accessor method for index file paths
  Returntype: string
  Exceptions: throws if file doesn't exist 
  Example   :

=cut



sub new_index{
  my ($self, $new_index) = @_;
  if($new_index){
    throw("SequenceIndexStatistics:new_index ".$new_index." should exist")
      unless(-e $new_index);
    $self->{new_index} = $new_index;
  }
  return $self->{new_index};
}

sub old_index{
  my ($self, $old_index) = @_;
  if($old_index){
    throw("SequenceIndexStatistics:old_index ".$old_index." should exist")
      unless(-e $old_index);
    $self->{old_index} = $old_index;
  }
  return $self->{old_index};
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


=head2 collection_type/name

  Arg [1]   : ReseqTrack::Tools::Statistics
  Arg [2]   : string, name or type of desired collection
  Function  : accessor method for collection name/type
  Returntype: string
  Exceptions: n/a
  Example   : n/a

=cut


sub collection_type{
  my ($self, $collection_type) = @_;
  if($collection_type){
    $self->{collection_type} = $collection_type;
  }
  return $self->{collection_type};
}

sub collection_name{
  my ($self, $collection_name) = @_;
  if($collection_name){
    $self->{collection_name} = $collection_name;
  }
  return $self->{collection_name};
}


=head2 run_id_hash

  Arg [1]   : ReseqTrack::Tools::Statistics
  Arg [2]   : hashref
  Function  : store run id hashref
  Returntype: hashref
  Exceptions: throws if not passed a hashref
  Example   :

=cut


sub run_id_hash{
  my ($self, $hash) = @_;
  if($hash){
    throw("Must pass run_id_hash a hashref not ".$hash) unless(ref($hash) eq 'HASH');
    $self->{run_id_hash} = $hash;
  }
  if(!$self->{run_id_hash}){
    $self->get_run_id_hash;
  }
  return $self->{run_id_hash};
}


=head2 get_run_id_list

  Arg [1]   : ReseqTrack::Tools::Statistics
  Arg [2]   : string, collection name
  Arg [3]   : string, colleciton type
  Arg [4]   : ReseqTrack::DBSQL::DBAdaptor;
  Function  : generate a hash of run ids on the basis of the given info
  Returntype: hashref
  Exceptions:
  Example   :

=cut



sub get_run_id_hash{
  my ($self, $collection_name, $collection_type, $db) = @_;
  my %hash;
  my $rmis = $self->rmis($collection_name, $collection_type, $db);
  foreach my $rmi(@$rmis){
    $hash{$rmi->run_id} = 1;
  }
  $self->{run_id_hash} = \%hash;
  return $self->run_id_hash;
}

sub rmis{
  my ($self, $rmis) = @_;
  if($rmis){
    $self->{rmis} = $rmis;
  }
  if(!$self->{rmis}){
    $self->get_rmis;
  }
  return $self->{rmis};
}

sub get_rmis{
  my ($self, $collection_name, $collection_type, $db) = @_;
  $collection_name = $self->collection_name unless($collection_name);
  $collection_type = $self->collection_type unless($collection_type);
  $db = $self->db unless($db);
   my @rmis;
  if($collection_name && $collection_type){
    my $ca = $db->get_CollectionAdaptor;
    my $collection =  $ca->fetch_by_name_and_type($collection_name, $collection_type);
    @rmis = @{$collection->others};
  }else{
    my $rmia = $db->get_RunMetaInfoAdaptor;
    @rmis = @{$rmia->fetch_all};
  }
  $self->{rmis} = \@rmis;
}


1;
