package ReseqTrack::DBSQL::BaseAdaptor;

require Exporter;
use vars qw(@ISA @EXPORT);
use strict;
use Benchmark;
use ReseqTrack::Tools::Exception qw(throw);
use DBI qw(:sql_types);

@ISA = qw(Exporter);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

sub new {
  my ($class,$dbobj) = @_;
  
  my $self = {};
  bless $self,$class;
  
  if( !defined $dbobj || !ref $dbobj ) {
    throw("Don't have a db [$dbobj] for new adaptor");
  }
  unless($dbobj->isa('ReseqTrack::DBSQL::DBAdaptor')){
    throw("Must be passed a ReseqTrack::DBSQL::DBAdaptor not a ".$dbobj);
  }
  $self->db($dbobj);
  $self->dbc($dbobj->dbc);
  return $self;
}

sub prepare{
  my ($self,$string) = @_;
 
# uncomment next line to cancel caching on the sql side. Needed for timing comparisons etc 
#  $string =~ s/SELECT/SELECT SQL_NO_CACHE/i;

  return $self->dbc->prepare($string);
}

=head2 db

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor $obj 
               the database this adaptor is using.
  Example    : $db = $adaptor->db();
  Description: Getter/Setter for the DatabaseConnection that this adaptor is 
               using.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
  Exceptions : none
  Caller     : Adaptors inherited fro BaseAdaptor
  Status     : Stable

=cut

sub db{
  my $self = shift;
  $self->{'db'} = shift if(@_);
  #print "Have db ".$self->{'db'}."\n";
  return $self->{'db'};

}

=head2 dbc

  Arg [1]    : (optional) Bio::EnsEMBL::DBSQL::DBConnection $obj 
               the database this adaptor is using.
  Example    : $db = $adaptor->db();
  Description: Getter/Setter for the DatabaseConnection that this adaptor is 
               using.
  Returntype : Bio::EnsEMBL::DBSQL::DBConnection
  Exceptions : none
  Caller     : Adaptors inherited fro BaseAdaptor
  Status     : Stable

=cut

sub dbc{
  my $self = shift;
  $self->{'dbc'} = shift if(@_);

  return $self->{'dbc'};
}

#Now standard fetch_all and fetch_by_dbID methods

sub columns{
  throw("Need to implement the columns method");
}

sub table_name{
  throw("Need to implement the table_name method");
}

sub where{
  return undef;
}

sub object_from_hashref{
  throw("Need to implement object_from_hashref method");
}

sub fetch_all{
  my ($self) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name;
  $sql .= " where ".$self->where if($self->where);
  my @objects;
  my $sth = $self->prepare($sql);
  eval{
    $sth->execute;
  };
  if($@){
    throw("Problem running $sql $@");
  }
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@objects, $object);
  }
  $sth->finish;
  return \@objects;
}


sub internal_id_column{
  my ($self) = @_;
  return $self->table_name."_id";
}

sub fetch_by_dbID{
   my ($self, $dbID) = @_;
   my $sql = "select ".$self->columns." from ".$self->table_name.
       " where ".$self->internal_id_column." = ?";   
   $sql .= " and ".$self->where if($self->where);
   my $sth = $self->prepare($sql);
   $sth->bind_param(1, $dbID);
   $sth->execute;
   my $rowHashref = $sth->fetchrow_hashref;
   my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
   $sth->finish;
   return $object;
}

sub fetch_by_column_name{
  my ($self, $column_name, $column_value) = @_;
  return $self->fetch_by_column_names([$column_name],[$column_value]);
}

sub fetch_by_column_names{
	my ($self, $column_names, $column_values) = @_;
	my @sql = ('select',$self->columns,' from ',$self->table_name, 'where 1 = 1');
	
	push @sql, "and", $self->where if ($self->where);
	
	for my $column_name (@$column_names){
		push @sql, 'and', $column_name, ' = ?';
	} 
	
	my $sth = $self->prepare(join ' ', @sql);
	
	my $index = 1;
	for my $value (@$column_values){
		$sth->bind_param($index++, $value);
	}
	
	$sth->execute;
	
	my @objects;
	while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
      push(@objects, $object);
    }
    $sth->finish;
    return \@objects;	
}

sub number_of_lines{
  my ($self) = @_;
  my $sql = "select count(*) from ".$self->table_name;
  my $sth = $self->prepare($sql);
  $sth->execute;
  my ($count) = $sth->fetchrow;
  $sth->finish;
  return $count;
}

sub store_history{
  my ($self, $object) = @_;
  throw("Can't store history for ".$object." that isnt a ReseqTrack::HasHistory")
      unless($object->isa("ReseqTrack::HasHistory"));
  my $hist_a = $self->db->get_HistoryAdaptor();
  if($object->history && @{$object->history} >= 1){
    foreach my $history(@{$object->history}){
      $history->other_id($object->dbID);
      $hist_a->store($history) unless($history->dbID);
    }
  }
}

sub store_attributes{
  my ($self, $object, $update) = @_;
  throw("Can't store statistics for ".$object." that isnt a ReseqTrack::HasHistory ")  unless($object->isa("ReseqTrack::HasHistory"));
  my $attr_a = $self->db->get_AttributeAdaptor();
  
  if($object->attributes){
    foreach my $statistics(@{$object->attributes}){
      $statistics->other_id($object->dbID);
      $statistics->table_name($object->object_table_name);
      $attr_a->store($statistics, $update);
    }
  }
}
1;
