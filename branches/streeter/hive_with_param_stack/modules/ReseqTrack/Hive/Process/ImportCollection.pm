
package ReseqTrack::Hive::Process::ImportCollection;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    $self->param_required('collection_type');
    $self->param_required('collection_name');
    $self->param_required('output_param');

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $ca = $db->get_CollectionAdaptor;
    my $collection = $ca->fetch_by_name_and_type($self->param('collection_name'), $self->param('collection_type'));

    if (!$collection) {
      $self->flows_non_factory(undef);
      return;
    }

    my $output_name_param = $self->param('output_param');
    my $output_values = [map {$_->name} @{$collection->others}];
    $self->output_param($output_name_param, $output_values);
}

1;

