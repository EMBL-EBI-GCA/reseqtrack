
package ReseqTrack::Hive::Process::SeedComplete;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $ps_id = $self->param_required('ps_id');

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $psa = $db->get_PipelineSeedAdaptor;
    my $pipeline_seed = $psa->fetch_by_dbID($ps_id);
    throw("did not get pipeline_seed for $ps_id") if !$pipeline_seed;

    # sanity check:
    my $self_url = $self->dbc->url;
    $self_url =~ s/[^\/]*\@//; # remove username and password from url
    my $ps_url = $pipeline_seed->hive_db->url;
    throw("url does not match $self_url $ps_url") if $self_url ne $ps_url;

    $psa->update_completed($pipeline_seed);
    
}

1;

