
package ReseqTrack::Hive::Process::ForeignFilesFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
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

    my $factory_timestamp = time();
    $self->output_param('factory_time', $factory_timestamp);

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

        my %output_params = (filename => $dropbox_files[0], dbmd5 => $file->md5,
                    dbsize => $file->size, dbID => $file->dbID);
        $self->prepare_factory_output( $file->dbID, \%output_params );
      }
    }
}

1;

