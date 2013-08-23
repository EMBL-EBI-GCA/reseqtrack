
package ReseqTrack::Hive::Process::SeedFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::GeneralUtils qw(delete_lock_string is_locked create_lock_string);
use Bio::EnsEMBL::Hive::Utils qw(destringify);

sub run {
    my $self = shift @_;

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $dbname = $self->dbc->dbname;
    my $host = $self->dbc->host;
    my $port = $self->dbc->port;
    my ($hive_db) = @{$db->get_HiveDBAdaptor->fetch_by_column_names(
            ['name', 'host', 'port'],
            [$dbname, $host, $port])};
    throw("did not get a hive_db object for $dbname $host $port") if !$hive_db;

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
    my %output_params;
    foreach my $param_name (@{$seeding_module->output_params}) {
      next if ! $self->param_is_defined($param_name);
      $output_params{$param_name} = $self->param($param_name);
    }

    my $psa = $db->get_PipelineSeedAdaptor;
    foreach my $seed_params (@{&$seeding_sub($pipeline, $seeding_options, \%output_params)}) {
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
      $self->prepare_factory_output_id($output_hash);
    }

    delete_lock_string($lock_string, $meta_adaptor);
    $hive_db->is_seeded(0);
    $db->get_HiveDBAdaptor->update($hive_db);
    
}

1;

