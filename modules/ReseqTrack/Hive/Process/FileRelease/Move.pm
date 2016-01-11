
package ReseqTrack::Hive::Process::FileRelease::Move;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists move_by_rsync);
use ReseqTrack::History;
use File::Copy qw(move);
use File::stat;
use File::Spec;
use File::Basename qw(fileparse);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub derive_path {
  my ($self, $dropbox_path, $file_object) = @_;
  my $derive_path_options = $self->param('derive_path_options');
  throw("Project-specific class must implement the derive_path subroutine");
}

sub param_defaults {
  return {
    'ps_attributes' => {},
    'derive_path_options' => {},
    'move_by_rsync' => 0,
    'collect' => 0,
  };
}

# An opportunity for the derive_path subroutine to make a last-second complaint.
sub reject_message {
  my ($self, $reject_message) = @_;
  if (defined $reject_message) {
    $self->param('reject_message', $reject_message);
  }
  return $self->param_is_defined('reject_message') ? $self->param('reject_message') : '';
}

sub run {
    my $self = shift @_;

    my $file_details = $self->param_required('file');
    my $db_params = $self->param_required('reseqtrack_db');
    my $hostname = $self->param_required('hostname');
    my $ps_attributes = $self->param('ps_attributes');
    my $collect = $self->param('collect');

    throw("input_id not correctly formulated") if ! defined $file_details->{'dropbox'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'db'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'dropbox'}->{'path'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'dropbox'}->{'ctime'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'db'}->{'dbID'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'db'}->{'updated'};

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});

    my $fa = $db->get_FileAdaptor;
    my $file_id = $file_details->{'db'}->{'dbID'};
    my $file_object = $fa->fetch_by_dbID($file_id);
    throw("did not find file $file_id in database") if !$file_object;

    my $index_object;
    my $collection_name = undef;
    my $index_extension = $file_details->{'dropbox'}->{'index_ext'};
    if ($index_extension) {
      my $index_file_id = $file_details->{'db'}->{'index_dbID'};
      $index_object = $fa->fetch_by_dbID($index_file_id);
      throw("did not find file $index_file_id in database") if !$index_object;
    }

    my $dropbox_path = $file_details->{'dropbox'}->{'path'};
    my $destination_path;
    ( $destination_path, $collection_name ) = $self->derive_path($dropbox_path, $file_object); 
    throw('collection name not found') if $collect && !$collection_name;
    
    if (my $reject_message = $self->reject_message) {
        $ps_attributes->{'message'} = $reject_message;
        $self->output_param('ps_attributes', $ps_attributes);
        $self->output_param('is_failed', 1);
        return;
    }

    if ($file_object->updated ne $file_details->{'db'}->{'updated'}
        || ($index_object && $index_object->updated ne $file_details->{'db'}->{'index_updated'})) {
        $ps_attributes->{'message'} = 'file updated in db since pipeline started';
        $self->output_param('ps_attributes', $ps_attributes);
        $self->output_param('is_failed', 1);
        return;
    }

    my $st = stat($dropbox_path) or throw("could not stat $dropbox_path: $!");
    my $index_st = $index_extension ? (stat($dropbox_path.$index_extension) or throw("could not stat $dropbox_path.$index_extension: $!"))
                : '';
    if ($st->ctime != $file_details->{'dropbox'}->{'ctime'}
        || ($index_extension && $index_st->ctime != $file_details->{'dropbox'}->{'index_ctime'})) {
        $ps_attributes->{'message'} = 'file changed since pipeline started';
        $self->output_param('ps_attributes', $ps_attributes);
        $self->output_param('is_failed', 1);
        return;
    }

    throw("derive_path subroutine did not return a full path: $destination_path")
        if !File::Spec->file_name_is_absolute($destination_path);
    my ($destination_filename, $destination_dir) = fileparse($destination_path);
    my ($dropbox_filename, $dropbox_dir) = fileparse($dropbox_path);
    
    my $change_name = $destination_filename ne $dropbox_filename ? 1 : 0;
    if ($change_name) {
      my $exists = $fa->fetch_by_filename($destination_filename);
      throw("file already exists in db: $dropbox_path $destination_filename" . $exists->[0]->name) if @$exists;
      if ($index_extension) {
        $exists = $fa->fetch_by_filename($destination_filename.$index_extension);
        throw("file already exists in db: $dropbox_path$index_extension $destination_filename$index_extension" . $exists->[0]->name) if @$exists;
      }
    }

    check_directory_exists($destination_dir);
    if ($self->param('move_by_rsync')) {
      $self->dbc->disconnect_when_inactive(1);
      $db->dbc->disconnect_when_inactive(1);
      move_by_rsync($dropbox_path, $destination_path);
      if ($index_extension) {
        move_by_rsync($dropbox_path.$index_extension, $destination_path.$index_extension);
      }
      $self->dbc->disconnect_when_inactive(0);
      $db->dbc->disconnect_when_inactive(0);
      throw("unexpected file size after rsync $dropbox_path $destination_path") if $file_object->size != -s $destination_path;
      throw("unexpected file size after rsync $dropbox_path $destination_path") if $index_extension && $index_object->size != -s $destination_path.$index_extension;
    }
    else {
      move($dropbox_path, $destination_path) or throw("error moving to $destination_path: $!");
      if( $index_extension ) {
        move($dropbox_path.$index_extension, $destination_path.$index_extension) or throw("error moving to $destination_path.$index_extension: $!");
      }
    }
    my $host = get_host_object($hostname, $db);
    $file_object->name($destination_path);
    $file_object->host($host);

    my $host_comment = 'changed host from '. $file_object->host->name.' to '. $host->name;
    my $host_history = ReseqTrack::History->new(
        -other_id => $file_id, -table_name => 'file', -comment => $host_comment);
    $file_object->history($host_history);
    if ($change_name) {
      my $name_comment = "changed filename from $dropbox_filename to $destination_filename";
      my $name_history = ReseqTrack::History->new(
          -other_id => $file_id, -table_name => 'file', -comment => $name_comment);
      $file_object->history($name_history);
      my $name_attribute = ReseqTrack::Attribute->new(
          -other_id => $file_id, -table_name => 'file',
          -attribute_name => 'ORIG_FILENAME',
          -attribute_value => $dropbox_filename);
      $file_object->attributes($name_attribute);
    }

    $fa->update($file_object,0,$change_name);
    
    if ($collect) {
       my $collection = ReseqTrack::Collection->new(
                          -name       => $collection_name, 
                          -type       => $file_object->type, 
                          -others     => $file_object, 
                          -table_name => 'file',
                        );

      $db->get_CollectionAdaptor->store($collection);
    }

    if ($index_extension) {
      $index_object->name($destination_path.$index_extension);
      $index_object->host($host);
      my $index_host_history = ReseqTrack::History->new(
          -other_id => $index_object->dbID, -table_name => 'file', -comment => $host_comment);
      $index_object->history($index_host_history);
      if ($change_name) {
        my $index_name_comment = "changed filename from $dropbox_filename$index_extension to $destination_filename$index_extension";
        my $index_name_history = ReseqTrack::History->new(
            -other_id => $file_id, -table_name => 'file', -comment => $index_name_comment);
        $file_object->history($index_name_history);
        my $index_name_attribute = ReseqTrack::Attribute->new(
            -other_id => $file_id, -table_name => 'file',
            -attribute_name => 'ORIG_FILENAME',
            -attribute_value => $dropbox_filename.$index_extension);
        $file_object->attributes($index_name_attribute);
      }
      $fa->update($index_object,0,$change_name);
    }
}

1;

