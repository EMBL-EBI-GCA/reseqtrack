
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

    my $dropbox_path = $file_details->{'dropbox'}->{'path'};
    my $destination_path = $self->derive_path($dropbox_path, $file_object);
    if (my $reject_message = $self->reject_message) {
        $ps_attributes->{'message'} = $reject_message;
        $self->output_param('ps_attributes', $ps_attributes);
        $self->output_param('is_failed', 1);
        return;
    }

    if ($file_object->updated ne $file_details->{'db'}->{'updated'}) {
        $ps_attributes->{'message'} = 'file updated in db since pipeline started';
        $self->output_param('ps_attributes', $ps_attributes);
        $self->output_param('is_failed', 1);
        return;
    }

    my $st = stat($dropbox_path) or throw("could not stat $dropbox_path: $!");
    if ($st->ctime != $file_details->{'dropbox'}->{'ctime'}) {
        $ps_attributes->{'message'} = 'file changed since pipeline started';
        $self->output_param('ps_attributes', $ps_attributes);
        $self->output_param('is_failed', 1);
        return;
    }

    throw("derive_path subroutine did not return a full path: $destination_path")
        if !File::Spec->file_name_is_absolute($destination_path);
    my ($destination_filename, $destination_dir) = fileparse($destination_path);
    if ($destination_filename ne fileparse($dropbox_path)) {
      my $exists = $fa->fetch_by_filename($destination_filename);
      throw("file already exists in db: " . $exists->[0]->name) if @$exists;
    }

    check_directory_exists($destination_dir);
    if ($self->param('move_by_rsync')) {
      $self->dbc->disconnect_when_inactive(1);
      $db->dbc->disconnect_when_inactive(1);
      move_by_rsync($dropbox_path, $destination_path);
      $self->dbc->disconnect_when_inactive(0);
      $db->dbc->disconnect_when_inactive(0);
      throw("unexpected file size after rsync $dropbox_path $destination_path") if $file_object->size != -s $destination_path;
    }
    else {
      move($dropbox_path, $destination_path) or throw("error moving to $destination_path: $!");
    }
    my $host = get_host_object($hostname, $db);
    my $comment = 'changed host from '. $file_object->host->name.' to '. $host->name;
    my $history = ReseqTrack::History->new(
        -other_id => $file_id, -table_name => 'file', -comment => $comment);
    $file_object->name($destination_path);
    $file_object->host($host);
    $file_object->history($history);
    $fa->update($file_object);
}

1;

