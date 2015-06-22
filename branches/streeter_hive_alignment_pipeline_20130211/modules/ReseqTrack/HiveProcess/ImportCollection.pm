
package ReseqTrack::HiveProcess::ImportCollection;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $collection_type = $self->param('collection_type') || die "'collection_type' is an obligatory parameter";
    my $collection_name = $self->param('collection_name') || die "'collection_name' is an obligatory parameter";

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $ca = $db->get_CollectionAdaptor;
    my $collection = $ca->fetch_by_name_and_type($collection_name, $collection_type);
    throw("Failed to find a collection for $collection_name $collection_type") if(!$collection);

    $self->output_this_branch('name' => [map {$_->name} @{$collection->others}]);
}

1;

