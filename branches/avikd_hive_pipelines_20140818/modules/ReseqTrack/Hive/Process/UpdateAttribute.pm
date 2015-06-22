package ReseqTrack::Hive::Process::UpdateAttribute;

use strict;
use warnings;
use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Collection;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::AttributeUtils;

sub param_defaults {
  return {
   attribute_prefix_column => 'CATEGORY',
  };
}

sub run {
    my $self = shift @_;
    
    $self->param_required('attribute_metrics');
    my $collection_type = $self->param_required('collection_type');
    my $collection_name = $self->param_required('collection_name');
    my $db_param = $self->param_required('reseqtrack_db');
    my $attribute_prefix_column = $self->param('attribute_prefix_column');
    
    my $attribute_metrics = $self->param_as_array('attribute_metrics'); 
 
    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_param});
    $db->dbc->disconnect_when_inactive(1);
    
    my $ca = $db->get_CollectionAdaptor;
    my $collection =  $ca->fetch_by_name_and_type( $collection_name, $collection_type );
    throw("Cannot find collection in db for $collection_name, $collection_type") if ( !$collection );
    
    my @statistics;
    
    for my $metrics_group (@$attribute_metrics) {
      $metrics_group = [$metrics_group] if( ref($metrics_group) eq 'HASH');
      warn ref($metrics_group),"\n";
      for my $metrics (@$metrics_group) {

        my $prefix = $metrics->{$attribute_prefix_column};
        while ( my ( $key, $value ) = each %$metrics ) {
          next if $key eq $attribute_prefix_column;
        
          if ($prefix) {
            $key = join( '_', $prefix, $key );
          }
          push @statistics, 
                 create_attribute_for_object( $collection, $key, $value )
                       if ( defined $value );
        }
      }
    }
    my $attributes = $collection->uniquify_attributes(\@statistics);
    $collection->attributes($attributes);
    $ca->store_attributes($collection);
}

1;
