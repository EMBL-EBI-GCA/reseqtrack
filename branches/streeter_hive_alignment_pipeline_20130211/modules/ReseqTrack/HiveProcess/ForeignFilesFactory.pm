
package ReseqTrack::HiveProcess::ForeignFilesFactory;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);
use File::Find qw(find);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});

    my $meta_adaptor = $db->get_MetaAdaptor;
    eval{ is_locked("file_release.lock", $meta_adaptor);};
    if ($@) {
      throw("ReseqTrack database is locked with file_release.lock in meta table");
    }
    create_lock_string("file_release.lock", $meta_adaptor);

    $self->output_this_branch(factory_timestamp => time());

    my $remote_hosts = $db->get_HostAdaptor->fetch_all_remote();
    my $fa = $db->get_FileAdaptor;
    foreach my $host (@$remote_hosts) {
      FILE:
      foreach my $file (@{$fa->fetch_by_host($host->dbID)}) {
        my $filename = $file->filename;
        my @dropbox_files;
        find( sub {push(@dropbox_files, $File::Find::name) if $_ eq $filename}, $host->dropbox_dir);
        throw("multiple files with the same name: @dropbox_files") if @dropbox_files >1;
        next FILE if !@dropbox_files;

        $self->output_child_branches(filename => $dropbox_files[0], md5 => $file->md5, size => $file->size, file_id => $file->dbID);
      }
    }
}

1;

