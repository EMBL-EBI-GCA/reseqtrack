
package ReseqTrack::Hive::PipeSeed::ForeignFiles;

use strict;
use warnings;
use base ('ReseqTrack::Hive::PipeSeed::BasePipeSeed');

use File::Find qw(find);
use File::stat;
use DateTime::Format::MySQL;

sub create_seed_params {
  my ($self) = @_;

  throw('this module will only accept pipelines that work on the file table')
      if $self->table_name ne 'file';

  my $db = $self->db;
  my $pipeline = $self->pipeline;

  my $remote_hosts = $db->get_HostAdaptor->fetch_all_remote();
  my $fa = $db->get_FileAdaptor;
  my $psa = $db->get_PipelineSeedAdaptor;
  my @seed_params;
  foreach my $host (@$remote_hosts) {
    FILE:
    foreach my $file (@{$fa->fetch_by_host($host->dbID)}) {
      my $existing_ps = $psa->fetch_by_seed_and_pipeline($file, $pipeline);
      if (@$existing_ps) {
        next FILE if grep {$_->is_running} @$existing_ps;
        next FILE if grep {$_->is_complete} @$existing_ps;
        next FILE if grep {$_->is_futile} @$existing_ps;
      }
      my $filename = $file->filename;
      my @dropbox_files;
      find( sub {push(@dropbox_files, $File::Find::name) if $_ eq $filename}, $host->dropbox_dir);
      throw("multiple files with the same name: @dropbox_files") if @dropbox_files >1;
      next FILE if !@dropbox_files;

      my $st = stat($dropbox_files[0]) or throw("could not stat $dropbox_files[0]: $!");
      my $drop_box_ctime = $st->ctime;
      my $updated = DateTime::Format::MySQL->parse_datetime($file->updated)->set_time_zone('local')->epoch;
      if ($drop_box_ctime > $updated) {
        $updated = $drop_box_ctime;
      }
      next FILE if grep {$_ > $updated} map {DateTime::Format::MySQL->parse_datetime($_->created)->set_time_zone('local')->epoch} @$existing_ps;

      my %dropbox_details = ('path' =>$dropbox_files[0], 'ctime' => $drop_box_ctime);
      my %db_details = (dbID => $file->dbID, updated => $file->updated);
      my %output_params = (file => {dropbox => \%dropbox_details, db => \%db_details});
      push(@seed_params, [$file, \%output_params]);
    }
  }
  $self->seed_params(\@seed_params);
}

1;
