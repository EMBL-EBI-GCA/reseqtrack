
package ReseqTrack::HiveProcess::FileReleaseMove;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use File::Copy qw(move);
use File::stat;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub derive_directory {
  throw("Project-specific class must implement the derive_directory subroutine");
}

sub run {
    my $self = shift @_;

    my $dropbox_filename = $self->param('dropbox_filename');
    my $hostname = $self->param('hostname') || '1000genomes.ebi.ac.uk';

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $fa = $db->get_FileAdaptor;

    my $file = $fa->fetch_by_dbID($self->param('db_file_id'));

    my $dir = $self->derive_directory;
    my $new_path = $dir . '/' . $file->filename;

    my $st = stat($dropbox_filename) or throw("could not stat $dropbox_filename: $!");
    if ($st->ctime > $self->param('branch_timestamp')) {
      my $log_adaptor = $db->get_RejectLogAdaptor;
      my $reject_log = ReseqTrack::RejectLog->new(
        -file_id => $self->param('db_file_id'),
        -is_reject => 'n',
        -reject_reason => "file changed since pipeline started",
        );
        $log_adaptor->store($reject_log, 1);
        return;
    }

    move($dropbox_filename, $new_path) or throw("error moving to $new_path: $!");
    my $host = get_host_object($hostname, $db);
    $file->name($new_path);
    $file->host($host);
    $fa->update($file);
}

1;

