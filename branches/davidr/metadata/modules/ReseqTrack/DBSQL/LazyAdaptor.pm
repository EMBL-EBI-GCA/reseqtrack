package ReseqTrack::DBSQL::LazyAdaptor;

use strict;
use base qw(ReseqTrack::DBSQL::BaseAdaptor);
use ReseqTrack::Tools::Exception qw(throw);

sub table_name{
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
		my $value = $hashref->{$column};
		
		
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
		if (defined $value) {
			push @bind_values,   $value;
			push @column_names,  $column;
			push @place_holders, '?';
		}
		elsif ($column eq 'created' || $column eq 'updated'){
			push @column_names, $column;
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
	
	my @bind_values;
	my @column_names;

	my $column_values = $self->column_values($obj);
	
	my $dbID_column = $self->internal_id_column;
	my $dbID = delete ${$column_values}{$dbID_column};
	
	
	while ( my ( $column, $value ) = each %$column_values ) {
		push @bind_values,   $value;
		push @column_names,  $column;
	}
	push @bind_values, $dbID;
	my $sql = "update ".$self->table_name." set ";
	
	while (my $col = shift @column_names){
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


sub columns {
	my ($self) = @_;
	my $tn = $self->table_name;
	return join ', ', map { "$tn.$_" } keys %{ $self->column_mappings( $self->object_class->new() ) };
}

1;
