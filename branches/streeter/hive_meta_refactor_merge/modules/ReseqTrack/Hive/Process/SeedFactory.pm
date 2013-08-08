
package ReseqTrack::Hive::Process::SeedFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::GeneralUtils qw(delete_lock_string is_locked create_lock_string);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $table_columns = $self->param_to_flat_array('table_column');
    my $filter_module = $self->param_is_defined('seed_filter_module') ? $self->param('seed_filter_module') : undef;
    my $filter_options = $self->param_is_defined('seed_filter_options') ? $self->param('seed_filter_options') : {};

    if (defined $filter_module) {
      eval "require $filter_module" or throw "cannot load module $module $@";
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

    my $table_name = $hive_db->pipeline->table_name;
    my $dbID_name = $table_name . '_id';

    my $adaptor = $db->get_adaptor_for_table($table_name);
    my $sql_existing =
          "SELECT $table_name.$dbID_name FROM $table_name, pipeline_seed, hive_db"
        . " WHERE pipeline_seed.hive_db_id = hive_db.hive_db_id"
        . " AND hive_db.pipeline_id = ?"
        . " AND (pipeline_seed.is_running = 1"
        .      " OR pipeline_seed.is_complete = 1"
        .      " OR pipeline_seed.is_futile = 1)";
    my $sql = "SELECT ".$adaptor->columns." FROM $table_name "
          . " WHERE $dbID_name NOT IN ($sql_existing)";
    if (my $type = $hive_db->pipeline->type) {
      $sql .= " AND $table_name.type = $type";
    }
    my $sth = $db->dbc->prepare($sql) or throw("could not prepare $sql: ".$db->dbc->errstr);
    $sth->bind_param(1, $hive_db->pipeline->dbID);
    $sth->execute or die "could not execute $sql: ".$sth->errstr;

    my @seed_ids;
    if (defined $filter_module) {
      my $filter_function = $filter_module . '::filter';
      SEED:
      while (my $rowHashref = $sth->fetchrow_hashref) {
        my $object = $adaptor->object_from_hashref($rowHashref);
        if (&$filter_function($object, $filter_options){
          push(@seed_ids, $object->dbID);
        }
      }
    }
    else {
      @seed_ids = map {$_->{$dbID_name}} @$fetchall_hashref;
    }

    my $psa = $db->get_PipelineSeedAdaptor;
    foreach my $seed_id (@seed_ids) {
      my $pipeline_seed = ReseqTrack::PipelineSeed->new
          (
           -seed_id         => $seed_id,
           -hive_db         => $hive_db,
           -is_running      => 1,
      );
      $psa->store($pipeline_seed);
      my %output_id = ('ps_id' => $pipeline_seed->dbID);
      foreach my $column_name (@$table_columns) {
        throw("$column_name is not a valid column name for $table_name") if !exists($rowHashref->{$column_name});
        $output_id{$column_name} = $rowHashref->{$column_name};
      }
      $self->prepare_factory_output_id($pipeline_seed->dbID, \%output_id);
    }
    $sth->finish;

    delete_lock_string($lock_string, $meta_adaptor);
    $hive_db->is_seeded(0);
    $db->get_HiveDBAdaptor->update($hive_db);
    
}

1;

