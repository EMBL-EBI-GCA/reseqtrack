
package ReseqTrack::Hive::Process::SeedFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::GeneralUtils qw(delete_lock_string is_locked create_lock_string);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub sql_existing {
  my ($self, $pipeline) = @_;
  my $table_name = $pipeline->table_name;
  my $dbID_name = $table_name . '_id';
  my $pipeline_id = $pipeline->dbID;
  my $sql_existing =
        "SELECT $table_name.$dbID_name FROM $table_name, pipeline_seed, hive_db"
      . " WHERE pipeline_seed.hive_db_id = hive_db.hive_db_id"
      . " AND hive_db.pipeline_id = $pipeline_id"
      . " AND (pipeline_seed.is_running = 1"
      .      " OR pipeline_seed.is_complete = 1"
      .      " OR pipeline_seed.is_futile = 1)";
  return $sql_existing;
}

sub run {
    my $self = shift @_;
    my $seed_labels = $self->param_to_flat_array('seed_label');

    if (!@$seed_labels) {
      $seed_labels = ['ps_id'];
    }

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $url = $self->dbc->url;
    $url =~ s/[^\/]*\@//; # remove username and password from url

    my $hive_db = $db->get_HiveDBAdaptor->fetch_by_url($url);
    throw("did not get a hive_db object for $url") if !$hive_db;

    my $meta_adaptor = $db->get_MetaAdaptor;
    my $lock_string = $hive_db->pipeline->name . '.lock';
    eval{ is_locked($lock_string, $meta_adaptor);};
    if ($@) {
      throw("ReseqTrack database is locked with $lock_string in meta table");
    }
    create_lock_string($lock_string, $meta_adaptor);

    my $pipeline = $hive_db->pipeline;
    my $seeding_module = $pipeline->seeding_module;
    eval "require $seeding_module" or throw "cannot load module $seeding_module $@";

    my $seeding_sub_name = $seeding_module . '::create_seed_params';
    my $seeding_sub = \&$seeding_sub_name;
    my $seeding_options = $pipeline->seeding_options 
                        ? destringify($pipeline->seeding_options)
                        : {};

    my $psa = $db->get_PipelineSeedAdaptor;
    foreach my $seed_params (@{&$seeding_sub($self, $pipeline, $seeding_options)}) {
      my ($seed, $output_hash) = @$seed_params;
      throw('seeding module has returned an object of the wrong type')
          if $seed->adaptor->table_name ne $pipeline->table_name;
      my $pipeline_seed = ReseqTrack::PipelineSeed->new
          (
           -seed_id         => $seed->dbID,
           -hive_db         => $hive_db,
           -is_running      => 1,
      );
      $psa->store($pipeline_seed);
      $output_hash->{'ps_id'} = $pipeline_seed->dbID;
      my @label_values;
      foreach my $seed_label (@$seed_labels) {
        throw("output hash is missing label $seed_label") if ! defined $output_hash->{$seed_label};
        push(@label_values, $output_hash->{$seed_label});
      }
      $self->prepare_factory_output_id(\@label_values, $output_hash);
    }

    delete_lock_string($lock_string, $meta_adaptor);
    $hive_db->is_seeded(0);
    $db->get_HiveDBAdaptor->update($hive_db);
    
}

1;

