
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

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $meta_adaptor = $db->get_MetaAdaptor;
    my $lock_string = 'pipeline.lock';
    eval{ is_locked($lock_string, $meta_adaptor);};
    if ($@) {
      throw("ReseqTrack database is locked with $lock_string in meta table");
    }
    create_lock_string($lock_string, $meta_adaptor);

    my $url = $self->dbc->url;

    my $hive_db = $db->get_HiveDBAdaptor->fetch_by_url($url);
    throw("did not get a hive_db object for $url") if !$hive_db;

    my $table_name = $hive_db->pipeline->table_name;
    my $table_column = $hive_db->pipeline->table_column;
    my $dbID_name = $table_name . '.' . $table_name . '_id';

    my $adaptor = $db->get_adaptor_for_table($table_name);
    my $sql_existing = "SELECT $dbID_name FROM $table_name, pipeline_seed, hive_db"
        . " WHERE pipeline_seed.hive_db_id = hive_db.hive_db_id"
        . " AND hive_db.pipeline_id = ?";
    my $sql = "SELECT ".$adaptor->columns." FROM $table_name "
          . " WHERE $dbID_name NOT IN ($sql_existing)";
    if (my $type = $hive_db->pipeline->type) {
      $sql .= " AND $table_name.type = $type";
    }
    my $sth = $db->dbc->prepare($sql) or throw("could not prepare $sql: ".$db->dbc->errstr);
    $sth->bind_param(1, $hive_db->pipeline->dbID);
    $sth->execute or die "could not execute $sql: ".$sth->errstr;;

    my $psa = $db->get_PipelineSeedAdaptor;
    while(my $rowHashref = $sth->fetchrow_hashref){
      my $seed_value = $rowHashref->{$table_name};
      my $pipeline_seed = ReseqTrack::PipelineSeed->new
          (
           -seed_id         =>$rowHashref->{$dbID_name},
           -hive_db         =>$hive_db
           -status          => 'RUNNING',
      );
      $psa->store($pipeline_seed);
      $self->prepare_factory_output_id($seed_value, {'ps_id' => $pipeline_seed->dbID, 'seed' => $seed_value});
    }
    $sth->finish;

    delete_lock_string($lock_string, $meta_adaptor);
    
}

1;

