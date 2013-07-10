
package ReseqTrack::Hive::Process::FileRelease::Move;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::History;
use File::Copy qw(move);
use File::stat;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub derive_directory {
  my ($dropbox_path, $file_object) = @_;
  throw("Project-specific class must implement the derive_directory subroutine");
}

sub run {
    my $self = shift @_;

    $self->param_required('file');
    my $file_details = $self->file_param('file');
    my $db_params = $self->param_required('reseqtrack_db');
    my $hostname = $self->param_required('hostname');

    throw("input_id not correctly formulated") if ! defined $file_details->{'dropbox'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'db'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'dropbox'}->{'path'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'dropbox'}->{'ctime'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'db'}->{'dbID'};
    throw("input_id not correctly formulated") if ! defined $file_details->{'db'}->{'updated'};

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});

    my $fa = $db->get_FileAdaptor;
    my $log_adaptor = $db->get_RejectLogAdaptor;
    my $file_id = $file_details->{'db'}->{'dbID'};
    my $file_object = $fa->fetch_by_dbID($file_id);
    throw("did not find file $file_id in database") if !$file_object;

    my $dropbox_path = $file_details->{'dropbox'}->{'path'};
    my $dir = $self->derive_directory($dropbox_path, $file_object);
    my $new_path = $dir . '/' . $file_object->filename;

    if ($file_object->updated ne $file_details->{'db'}->{'updated'}) {
      my $reject_log = ReseqTrack::RejectLog->new(
        -file_id => $file_id,
        -is_reject => 'n',
        -reject_reason => "job is out of date with database",
        );
        $log_adaptor->store($reject_log, 1);
        $self->flows_non_factory(undef);
        return;
    }

    my $st = stat($dropbox_path) or throw("could not stat $dropbox_path: $!");
    if ($st->ctime != $file_details->{'dropbox'}->{'ctime'}) {
      my $reject_log = ReseqTrack::RejectLog->new(
        -file_id => $file_id,
        -is_reject => 'n',
        -reject_reason => "file changed since pipeline started",
        );
        $log_adaptor->store($reject_log, 1);
        $self->flows_non_factory(undef);
        return;
    }

    move($dropbox_path, $new_path) or throw("error moving to $new_path: $!");
    my $host = get_host_object($hostname, $db);
    my $comment = 'changed host from '. $file_object->host->name.' to '. $host->name;
    my $history = ReseqTrack::History->new(
        -other_id => $file_id, -table_name => 'file', -comment => $comment);
    $file_object->name($new_path);
    $file_object->host($host);
    $file_object->history($history);
    $fa->update($file_object);
    $log_adaptor->clear_table_by_file_id($file_id);
}

1;

