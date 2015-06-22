=pod

=head1 NAME

ReseqTrack::Tools::CollectionUtils

=head1 SYNOPSIS

This is a collection of methods to standardise Collection object handling

=head1 Example

use ReseqTrack::Tools::CollectionUtils qw(copy_collection_object);

my $new_collection = copy_collection_object($new_collection);

=cut

package ReseqTrack::Tools::CollectionUtils;

use strict;
use Exporter;

use strict;
use Exporter;
use ReseqTrack::Collection;
use ReseqTrack::History;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(copy_collection_object calculate_collection_comment create_collection_history_object);


=head2 copy_collection_object

  Arg [1]   : ReseqTrack::Collection
  Function  : to create a copy of the give collection object
  Returntype: ReseqTrack::Collection
  Exceptions: 
  Example   : 

=cut


sub copy_collection_object{
  my ($collection) = @_;
  my $new_collection = ReseqTrack::Collection->new(
    -dbID => $collection->dbID,
    -type => $collection->type,
    -name => $collection->name,
    -others => $collection->others,
    -table_name => $collection->table_name,
      );
  return $new_collection;
}


=head2 create_collection_history_object

  Arg [1]   : ReseqTrack::Collection
  Arg [2]   : ReseqTrack::Collection
  Arg [3]   : string, comment
  Function  : create history object, if necessary work out the difference to comment
  on
  Returntype: ReseqTrack::History
  Exceptions: throws if old object (second arg) doesn't have a defined dbID
  Example   : my $history = create_collection_history(undef, $col, "CHANGED TYPE");

=cut


sub create_collection_history_object{
  my ($new, $old, $comment) = @_;
  throw("ReseqTrack::Tools::CollectionUtils ".$old." needs a dbID otherwise you ".
        "can't create a history object for it") unless($old->dbID);
  $comment = calculate_collection_comment($new, $old) unless($comment);
  if($comment){
    my $history = ReseqTrack::History->new(
      -other_id => $old->dbID,
      -table_name => 'collection',
      -comment => $comment,
        );
    return $history;
  }else{
    my ($p, $f, $l) = caller;
    print STDERR "Can't find a difference between ".$new." and ".$old." $p $l\n";
    return undef;
  }

}


=head2 calculate_collection_comment

  Arg [1]   : ReseqTrack::Collection
  Arg [2]   : ReseqTrack::Collection
  Function  : finds the difference between the two objects
  Returntype: string, comment info
  Exceptions: throws if name or table name are different
  Example   : $comment = calculate_collection_comment($new, $old);

=cut


sub calculate_collection_comment{
  my ($new, $old) = @_;
  throw("Need both an old and a new object to calculate a comment about")
      unless($new && $old);
  if($new->name ne $old->name){
    throw($old." ".$old->name." does not match ".$new." ".$new->name." are you ".
          "comparing the correct objects");
  }
  if($new->table_name ne $old->table_name){
    throw("You can't change the table associated with ".$old->name." from ".
          $old->table_name." to ".$new->table_name);
  }
  if($new->type ne $old->type){
    return "Changing collection type from ".$old->type." to ".$new->type;
  }
  return undef;
}

1;

