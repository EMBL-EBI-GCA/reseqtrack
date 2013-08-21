
package ReseqTrack::Hive::Process::UpdateSeed;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Attribute;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub param_defaults {
  return {
    is_complete => 0,
    is_failed => 0,
    is_futile => 0,
    delete_seeds => 0,
    ps_attributes => {}
  };
}

sub run {
    my $self = shift @_;

    my $ps_id = $self->param_required('ps_id');
    my $is_complete = $self->param('is_complete') ? 1 : 0;
    my $is_failed = $self->param('is_failed') ? 1 : 0;
    my $is_futile = $self->param('is_futile') ? 1 : 0;
    my $delete_seeds = $self->param('delete_seeds') ? 1 : 0;
    my $attributes = $self->param('ps_attributes');
    $is_failed ||= $is_futile;
    throw('one of the following must be set: is_failed, is_complete, delete_seeds')
        if !grep {$_} ($is_failed, $is_complete, $delete_seeds);

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $psa = $db->get_PipelineSeedAdaptor;
    my $pipeline_seed = $psa->fetch_by_dbID($ps_id);
    throw("did not get pipeline_seed for $ps_id") if !$pipeline_seed;

    # sanity check:
    my $self_dbname = $self->dbc->dbname;
    my $ps_dbname = $pipeline_seed->hive_db->name;
    throw("dbnames do not match $self_dbname $ps_dbname") if $self_dbname ne $ps_dbname;

    if ($delete_seeds) {
      my $all_seeds = $psa->fetch_by_seed_and_pipeline($pipeline_seed->seed, $pipeline_seed->pipeline);
      $psa->delete($all_seeds);
      return;
    }

    while (my ($key, $value) = each %$attributes) {
      my $attribute = ReseqTrack::Attribute->new(
        -table_name => 'pipeline_seed', -other_id => $ps_id,
        -attribute_name => $key,
        -attribute_value => $value,
        );
      $pipeline_seed->attributes($attribute);
    }

    if ($is_failed) {
      $psa->update_failed($pipeline_seed, $is_futile);
    }
    else {
      $psa->update_completed($pipeline_seed);
    }
    
}

1;

