
package ReseqTrack::Hive::Process::UpdateSeed;

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
    my $is_complete = ($self->param_is_defined('is_complete') && $self->param('is_complete')) ? 1 : 0;
    my $is_failed = ($self->param_is_defined('is_failed') && $self->param('is_failed')) ? 1 : 0;
    my $is_futile = ($self->param_is_defined('is_futile') && $self->param('is_futile')) ? 1 : 0;
    $is_failed ||= $is_futile;
    throw('is_complete or is_failed must be set') if !$is_failed && !$is_complete;

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $psa = $db->get_PipelineSeedAdaptor;
    my $pipeline_seed = $psa->fetch_by_dbID($ps_id);
    throw("did not get pipeline_seed for $ps_id") if !$pipeline_seed;

    # sanity check:
    my $self_url = $self->dbc->url;
    $self_url =~ s/[^\/]*\@//; # remove username and password from url
    my $ps_url = $pipeline_seed->hive_db->url;
    throw("url does not match $self_url $ps_url") if $self_url ne $ps_url;

    if ($is_failed) {
      $psa->update_failed($pipeline_seed, $is_futile);
    }
    else {
      $psa->update_completed($pipeline_seed);
    }
    
}

1;

