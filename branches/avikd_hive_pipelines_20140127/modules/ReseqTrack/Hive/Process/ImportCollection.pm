
package ReseqTrack::Hive::Process::ImportCollection;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);


sub param_defaults {
  return {
    collection_type => undef,
    collection_name => undef,
    collection_id => undef,
  };
}

sub run {
    my $self = shift @_;

    my $collection_type = $self->param('collection_type');
    my $collection_name = $self->param('collection_name');
    my $collection_id = $self->param('collection_id');
    $self->param_required('output_param');

    throw("must have a collection name or collection id")
          if (!defined $collection_id && !defined $collection_name);
    throw("must have a collection type with a collection name")
          if (defined $collection_name && !defined $collection_type);

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param_required('reseqtrack_db')});
    my $ca = $db->get_CollectionAdaptor;
    my $collection = $ca->fetch_by_name_and_type($self->param('collection_name'), $self->param('collection_type'));

    my $collection = defined $collection_id ? $ca->fetch_by_dbID($collection_id)
                  : $ca->fetch_by_name_and_type($collection_name, $collection_type);

    my $output_name_param = $self->param('output_param');
    my $output_values = $collection ? [map {$_->name} @{$collection->others}] : [];
    $self->output_param($output_name_param, $output_values);
}

1;

