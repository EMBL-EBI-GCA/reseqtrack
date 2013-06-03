
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

    $self->param_required('collection_type');
    $self->param_required('collection_name_param');
    $self->param_required('output_name_param');

    my $param_name = $self->param('collection_name_param');
    $self->param_required($param_name);

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $ca = $db->get_CollectionAdaptor;
    my $collection = $ca->fetch_by_name_and_type($self->param($param_name), $self->param('collection_type')));
    throw(join(' ', 'Failed to find a collection for',$self->param($param_name), $self->param('collection_type')))if(!$collection);

    my $output_name_param = $self->param('output_name_param');
    my $output_values = map {$_->name} @{$collection->others};
    if ($collection->table_name eq 'file') {
      $output_values = $self->make_file_ids($output_values);
    }
    $self->param($output_name_param, $output_values);
    $self->add_to_dataflow($output_name_param);
}

1;

