
package ReseqTrack::Hive::Process::LoadFile;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::FileUtils qw(create_object_from_path assign_type assign_type_by_filename);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Collection;
use ReseqTrack::History;
use Digest::MD5::File qw(file_md5_hex);


sub param_defaults {
  return {
    type => undef,
    md5 => {},
    collect => 0,
    host_name => '1000genomes.ebi.ac.uk',
  };
}


sub run {
    my $self = shift @_;

    $self->param_required('file');
    my $db_params = $self->param_required('reseqtrack_db');
    my $type = $self->param('type');
    my $md5_hash = $self->file_param('md5');
    my $collect = $self->param('collect') ? 1 : 0;
    my $collection_name = $collect ? $self->param_required('collection_name') : undef;
    my $host_name = $self->param('host_name');
    my $pipeline_seed_id = $self->param_required('ps_id');

    my $file_paths = $self->file_param_to_flat_array('file');
    throw('no files') if !@$file_paths;

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});
    my $fa = $db->get_FileAdaptor;
    my $poa = $db->get_PipelineOutputAdaptor;
    my $host = get_host_object($host_name, $db);

    my $pipeline_seed;

    my @files;
    my $file_type_rules = defined $type ? undef : $db->get_FileTypeRuleAdaptor->fetch_all_in_order;
    foreach my $path (@$file_paths) {
      check_file_exists($path);
      $type //= assign_type_by_filename($path, $file_type_rules);
      throw("could not assign a type for $path") if !defined $type || $type eq '';
      my $file = create_object_from_path($path, $type, $host);
      my $md5 = $md5_hash->{$file};
      if (! defined $md5) {
        $db->dbc->disconnect_when_inactive(1);
        $md5 = file_md5_hex($path);
        $db->dbc->disconnect_when_inactive(0);
      }
      $file->md5($md5);
      
      my $existing_file = $fa->fetch_by_name($path);
      if ($existing_file) {
        $pipeline_seed //= $db->get_PipelineSeedAdaptor->fetch_by_dbID($pipeline_seed_id);
        $file->dbID($existing_file->dbID);
        my $history = ReseqTrack::History->new(
          -other_id => $file->dbID, -table_name => 'file',
          -comment => 'updated by pipeline '.$pipeline_seed->pipeline->name,
        );
        $file->history($history);
        $fa->update($file);
      }
      else {
        $fa->store($file);
      }

      my $pipeline_output = ReseqTrack::PipelineOutput->new(
        -pipeline_seed_id => $pipeline_seed_id,
        -table_name => 'file', -output => $file,
        -action => 'CREATED',
      );
      $poa->store($pipeline_output);
      push(@files, $file);
    }

    if ($collect) {
      my $collection = ReseqTrack::Collection->new(
          -name => $collection_name, -type => $type // $files[0]->type,
          -others => \@files, -table_name =>'file',
      );
      $db->get_CollectionAdaptor->store($collection);
    }
}

1;

