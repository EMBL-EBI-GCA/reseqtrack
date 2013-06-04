package ReseqTrack::DBSQL::StudyAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::LazyAdaptor);
use ReseqTrack::Study;

sub new {
	my ( $class, $db ) = @_;

	my $self = $class->SUPER::new($db);

	return $self;
}

sub column_mappings {
	my ($self,$s) = @_;
	
	throw("must be passed an object") unless ($s);
	
	return {
		study_id     => sub { $s->dbID(@_) } ,
		status       => sub { $s->status(@_) },
		md5          => sub { $s->md5(@_) },
		type         => sub { $s->type(@_) },
		title        => sub { $s->title(@_) },
		source_id		 => sub { $s->source_id(@_) },
	};
}

sub object_class {
	return 'ReseqTrack::Study';
}

sub table_name {
	return "study";
}

sub store {
	my ( $self, $study, $update ) = @_;
	my $existing_record = $self->fetch_by_dbID( $study->dbID ) if ($study->dbID);

	if ( $existing_record ) {
		if ($update){
			return $self->update($study);	
		}
	}
	
	$self->SUPER::store($study,$update);

}

sub fetch_by_source_id{
   my ($self, $source_id) = @_;
   return pop @{$self->fetch_by_column_name("source_id",$source_id)};
}


1;
