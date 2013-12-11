
package ReseqTrack::Hive::Process::LoadFile;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::FileUtils qw(create_object_from_path assign_type assign_type_by_filename);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_directory_exists check_file_does_not_exist);
use ReseqTrack::Collection;
use ReseqTrack::History;
use Digest::MD5::File qw(file_md5_hex);
use File::Copy qw( move );
use File::Basename qw( dirname );
use Cwd qw( abs_path);


sub param_defaults {
  return {
    type => undef,
    md5 => {},
    collect => 0,
    name_file_module => 'ReseqTrack::Hive::NameFile::BaseNameFile',
    name_file_method => 'default',
    name_file_params => {},
    clobber => 0,
    host_name => '1000genomes.ebi.ac.uk',
  };
}


sub run {
    my $self = shift @_;

    $self->param_required('file');
    my $db_params = $self->param_required('reseqtrack_db');
    my $type = $self->param('type');
    my $md5_hash = $self->param('md5');
    my $collect = $self->param('collect') ? 1 : 0;
    my $collection_name = $collect ? $self->param_required('collection_name') : undef;
    my $host_name = $self->param('host_name');
    my $pipeline_seed_id = $self->param_required('ps_id');
    my $clobber = $self->param('clobber') ? 1 : 0;

    my $current_file_paths = $self->param_as_array('file');
    throw('no files') if !@$current_file_paths;

    foreach my $current_path (@$current_file_paths) {
      check_file_exists($current_path);
    }

    my $name_file_module = $self->param('name_file_module') // param_defaults()->{'name_file_module'};
    eval "require $name_file_module" or throw "cannot load module $name_file_module $@";
    my $renamer = $name_file_module->new(
          -method => $self->param('name_file_method') // param_defaults()->{'name_file_method'},
          -params => $self->param('name_file_params'),
          -files => $current_file_paths,
          );
    $renamer->derive_paths;
    my $file_map = $renamer->file_map;

    my $load_file_paths = [];
    FILE:
    foreach my $old_path (@$current_file_paths) {
      my $new_path = $file_map->{$old_path};
      if($new_path) {;
        my $dir_name = dirname($new_path);
        check_directory_exists($dir_name);
        if (!$clobber && (abs_path($new_path) ne abs_path($old_path))) {
          check_file_does_not_exist($new_path);
        }
      }
      push(@$load_file_paths, $new_path || $old_path);
    }

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});

    my %types;
    if (!defined $type) {
      my $file_type_rules = $db->get_FileTypeRuleAdaptor->fetch_all_in_order;
      foreach my $load_path (@$load_file_paths) {
        my $assigned_type = assign_type_by_filename($load_path, $file_type_rules);
        throw("could not assign a type for $load_path") if !defined $assigned_type || $assigned_type eq '';
        $types{$load_path} = $assigned_type;
      }
    }

    if (grep {!$md5_hash->{$_}} @$current_file_paths) {
      $db->dbc->disconnect_when_inactive(1);
      foreach my $current_path (@$current_file_paths) {
        $md5_hash->{$current_path} //= file_md5_hex($current_path);
      }
      $db->dbc->disconnect_when_inactive(0);
    }

    my $fa = $db->get_FileAdaptor;
    my $poa = $db->get_PipelineOutputAdaptor;
    my $host = get_host_object($host_name, $db);
    my $pipeline_seed;


    my @file_objects;
    foreach my $current_path (@$current_file_paths) {
      my $load_path = $file_map->{$current_path} || $current_path;
      my $file = create_object_from_path($load_path, $type // $types{$load_path}, $host);
      my $md5 = $md5_hash->{$current_path};
      $file->md5($md5);

      if (my $move_path = $file_map->{$current_path}) {
        move($current_path, $move_path) or throw("could not move $current_path to $move_path $!");
      }
      
      my $action;
      if (my $existing_file = $fa->fetch_by_name($load_path)) {
        $pipeline_seed //= $db->get_PipelineSeedAdaptor->fetch_by_dbID($pipeline_seed_id);
        $file->dbID($existing_file->dbID);
        my $history = ReseqTrack::History->new(
          -other_id => $file->dbID, -table_name => 'file',
          -comment => 'updated by pipeline '.$pipeline_seed->pipeline->name,
        );
        $file->history($history);
        $fa->update($file);
        $action = 'UPDATED';
      }
      else {
        $fa->store($file);
        $action = 'CREATED';
      }

      my $pipeline_output = ReseqTrack::PipelineOutput->new(
        -pipeline_seed_id => $pipeline_seed_id,
        -table_name => 'file', -output => $file,
        -action => $action,
      );
      $poa->store($pipeline_output);
      push(@file_objects, $file);
    }

    if ($collect) {
      my $collection = ReseqTrack::Collection->new(
          -name => $collection_name, -type => $type // $file_objects[0]->type,
          -others => \@file_objects, -table_name =>'file',
      );
      $db->get_CollectionAdaptor->store($collection);
    }

    if (@$load_file_paths == 1) {
      $self->output_param('file', $load_file_paths->[0]);
    }
    else {
      $self->output_param('file', $load_file_paths);
    }

}

1;

