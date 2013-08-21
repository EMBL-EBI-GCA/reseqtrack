package ReseqTrack::DBSQL::LazyAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::BaseAdaptor);
use ReseqTrack::Tools::Exception qw(throw);
use List::Compare;

sub table_name {
  throw("Need to implement the table_name method");
}

=head2
 Arg[1] : Object of the class handled by this adaptor
 Returntype: hashref. keys are column names, values are sub routines for getter/setter methods.
 Example implementation:
 	my ($self,$s) = @_;
	
	throw("must be passed an object") unless ($s);
	
	return {
		study_id     => sub { $s->dbID(@_) } ,
		ena_study_id => sub { $s->ena_study_id(@_) },
		status       => sub { $s->status(@_) },
		md5          => sub { $s->md5(@_) },
		type         => sub { $s->type(@_) },
		title        => sub { $s->title(@_) },
	};

=cut

sub column_mappings {
  throw("Need to implement the column_mappings method");
}

sub object_class {
  throw("Need to specify the class this adaptor returns ");
}

sub object_from_hashref {
  my ( $self, $hashref ) = @_;
  throw("Can't create an object from an empty hashref")
    unless ($hashref);

  my $obj = $self->object_class->new;
  $obj->adaptor($self);

  my $column_mappings = $self->column_mappings($obj);

  for my $column ( keys %$column_mappings ) {
    my $method = $column_mappings->{$column};
    my $value  = $hashref->{$column};

    $method->($value);
  }

  return $obj;
}

sub column_values {
  my ( $self, $obj ) = @_;
  throw("must be passed an object") unless ($obj);
  my %col_val;

  my $column_mappings = $self->column_mappings($obj);

  for my $column ( keys %$column_mappings ) {
    my $method = $column_mappings->{$column};
    my $value  = $method->();

    $col_val{$column} = $value;
  }
  return \%col_val;
}

=head2 store
	Inserts the passed object to the db.
	Returns the object with the dbID set.
=cut

sub store {
  my ( $self, $obj ) = @_;
  throw("must be passed an object") unless ($obj);
  my @bind_values;
  my @column_names;
  my @place_holders;

  my $column_values = $self->column_values($obj);

  while ( my ( $column, $value ) = each %$column_values ) {
    if ( defined $value ) {
      push @bind_values,   $value;
      push @column_names,  $column;
      push @place_holders, '?';
    }
    elsif ( $column eq 'created' || $column eq 'updated' ) {
      push @column_names,  $column;
      push @place_holders, 'now()';
    }
  }

  my $sql =
      "insert into "
    . $self->table_name . " ("
    . join( ",", @column_names )
    . ") values ("
    . join( ",", @place_holders ) . ")";

  my $sth = $self->prepare($sql);

  my $param_index = 0;
  for my $val (@bind_values) {
    $sth->bind_param( ++$param_index, $val );
  }

  $sth->execute;

  my $dbID = $sth->{'mysql_insertid'};

  $sth->finish;
  $obj->dbID($dbID);
  $obj->adaptor($self);
  return $obj;
}

sub update {
  my ( $self, $obj ) = @_;

  $self->write_history($obj);

  my @bind_values;
  my @column_names;

  my $column_values = $self->column_values($obj);

  my $dbID_column = $self->internal_id_column;
  my $dbID        = delete ${$column_values}{$dbID_column};

  while ( my ( $column, $value ) = each %$column_values ) {
    push @bind_values,  $value;
    push @column_names, $column;
  }
  push @bind_values, $dbID;
  my $sql = "update " . $self->table_name . " set ";

  while ( my $col = shift @column_names ) {
    $sql .= "$col = ?";
    $sql .= ', ' if (@column_names);
  }
  $sql .= " where $dbID_column = ?";

  my $sth = $self->prepare($sql);

  my $param_index = 0;
  for my $val (@bind_values) {
    $sth->bind_param( ++$param_index, $val );
  }

  $sth->execute;
  $sth->finish;
  return $obj;
}

sub find_column_changes {
  my ($self,$new_obj,$old_obj) = @_;
  
  my @changes;
  my $old_values = $self->column_values($old_obj);
  my $new_values = $self->column_values($new_obj);

  for my $field ( keys %$new_values ) {
    my $ov = ( $old_values->{$field} || '<null>' );
    my $nv = ( $new_values->{$field} || '<null>' );

    if ( $ov ne $nv ) {
      push @changes, "$field:$ov to $nv";
    }
  }
  return @changes;
}

sub find_attribute_changes {
  my ($self,$new_obj,$old_obj) = @_;
  
  my @changes;
  
  my $oah = $old_obj->attributes_hash;
  my $nah = $new_obj->attributes_hash;

  my @old_keys = sort keys %$oah;
  my @new_keys = sort keys %$nah;

  my $lc      = List::Compare->new( \@old_keys, \@new_keys );
  my @removed = $lc->get_Lonly;
  my @added   = $lc->get_Ronly;

  if ( $lc->get_Lonly ) {
    my $comment = '-attrs:';
    $comment .= join( ',', $lc->get_Lonly );
    push @changes, $comment;
  }
  if ( $lc->get_Ronly ) {
    my $comment = '+attrs:';
    $comment .= join( ',', $lc->get_Ronly );
    push @changes, $comment;
  }
  
  for my $field ( $lc->get_intersection ) {
    my $ov = ( $oah->{$field}->attribute_value || '<null>' );
    my $nv = ( $nah->{$field}->attribute_value || '<null>' );

    if ( $ov ne $nv ) {
      my $comment = join( ' ', 'attr', $field, $ov, 'to', $nv );
      push @changes, $comment;
    }
  }
  
  return @changes;
}

sub write_history {
  my ( $self, $new_obj ) = @_;
  my $old_obj = $self->fetch_by_dbID( $new_obj->dbID );

  my @changes;
  push @changes, $self->find_column_changes($new_obj,$old_obj);
  push @changes, $self->find_attribute_changes($new_obj,$old_obj);
  
  if (@changes) {
    my $history = ReseqTrack::History->new(
      -other_id   => $new_obj->dbID,
      -table_name => $self->table_name,
      -comment    => join( ';', @changes ),
    );
    $new_obj->history($history);
    $self->store_history($new_obj);
  }
}

sub columns {
  my ($self) = @_;
  my $tn = $self->table_name;
  return join ', ', map { "$tn.$_" }
    keys %{ $self->column_mappings( $self->object_class->new() ) };
}

1;
