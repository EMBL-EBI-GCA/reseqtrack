package ReseqTrack::DBSQL::ERA::StudyAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::ERA::BaseEraAdaptor);
use ReseqTrack::Tools::Exception qw(throw);

sub column_mappings {
	my ($self,$s) = @_;
	
	throw("must be passed an object") unless ($s);
	
	return {
		study_id     => sub { $s->source_id(@_) } ,
		status       => sub { $s->status(@_) },
		md5          => sub { $s->md5(@_) },
		type         => sub { $s->type(@_) },
		title        => sub { $s->title(@_) },
	};
}

sub object_class {
	return 'ReseqTrack::Study';
}

sub table_name {
	return "study";
}

sub fetch_by_source_id{
   my ($self, $source_id) = @_;
   return pop @{$self->fetch_by_column_name("study_id",$source_id)};
}

sub xml_column {
	return "";
}

1;