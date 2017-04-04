
package ReseqTrack::Hive::Process::FileRelease::RemoveLock;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::GeneralUtils qw(delete_lock_string);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $meta_adaptor = $db->get_MetaAdaptor;
    delete_lock_string("file_release.lock", $meta_adaptor);
    $db->dbc->disconnect_if_idle();
}

1;

