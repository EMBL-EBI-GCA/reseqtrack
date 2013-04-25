
package ReseqTrack::HiveProcess::FileReleaseMove;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use File::Copy qw(move);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub derive_directory {
  throw("Project-specific class must implement the derive_directory subroutine");
}

sub run {
    my $self = shift @_;

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $fa = $db->get_FileAdaptor;

    my $file = $fa->fetch_by_dbID($self->param('db_file_id'));

    my $dir = $self->derive_directory;
    my $new_path = $dir . '/' . $file->filename;
    move($self->param('drobox_filename'), $new_path) or throw("error moving to $new_path: $!");
    $file->name($new_path);
    $fa->update($file);
}

1;

