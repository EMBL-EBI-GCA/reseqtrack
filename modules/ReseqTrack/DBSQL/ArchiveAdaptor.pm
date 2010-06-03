package ReseqTrack::DBSQL::ArchiveAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Archive;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::GeneralUtils;

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db, $nolock) = @_;

  my $self = $class->SUPER::new($db);
  my $ma = $self->db->get_MetaAdaptor;
  unless($nolock){
    unless(is_locked("archive.lock", $ma)){
      create_lock_string("archive.lock", $ma);
    }
  }
  return $self;
}

sub columns{
  my ($self) = @_;
  return "archive.archive_id, archive.name, archive.file_id, archive.md5, archive.size, archive.relative_path, archive.volume_name, archive.created, archive.updated, archive.priority, archive.new_name, archive.new_relative_path, archive.archive_action_id, archive.archive_location_id, archive.fire_action_id, archive.fire_exit_code, archive.fire_exit_reason";
}

sub table_name{
  return "archive";
}


sub fetch_by_name{
  my ($self, $name) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where ".$self->table_name.".name = ?";
  my $sth = $self->prepare;
  $sth->bind_param(1, $name);
  $sth->execute;
  my @archives;
  while(my ($rowHashref) = $sth->fetchrow_hashref){
    push(@archives, $self->object_from_hashref($rowHashref));
  }
  if(@archives >= 2){
    throw("Can't deal with two or more archive objects with the same name");
  }
  return $archives[0];
}


sub store{
  my ($self, $archive) = @_;
  throw("Can't store ".$archive." using ReseqTrack::DBSQL::ArchiveAdaptor") 
      unless($archive->isa("ReseqTrack::Archive"));
  throw("Can't archive ".$archive->file->full_path." which doesn't have an md5 defined")
      unless($archive->file->md5);
  throw("Can't store archive without an archive_action_id and an archive_location id") unless($archive->archive_action_id && $archive->archive_location_id);
  my $sql = "insert ignore into archive (name, file_id, md5, size, relative_path, ".
      "volume_name, created, updated, priority, new_name, new_relative_path, archive_action_id, archive_location_id) ".
      "values(?, ?, ?, ?, ?, ?, now(), now(), ?, ?, ?, ?, ?)";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $archive->name);
  $sth->bind_param(2, $archive->file_id);
  $sth->bind_param(3, $archive->md5);
  $sth->bind_param(4, $archive->size);
  $sth->bind_param(5, $archive->relative_path);
  $sth->bind_param(6, $archive->volume_name);
  $sth->bind_param(7, $archive->priority);
  $sth->bind_param(8, $archive->new_name);
  $sth->bind_param(9, $archive->new_relative_path);
  $sth->bind_param(10, $archive->archive_action_id);
  $sth->bind_param(11, $archive->archive_location_id);
  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish();
  $archive->dbID($dbID),
  return $archive;
}



sub remove{
  my ($self, $archive) = @_;
   throw("Can't update ".$archive." using ReseqTrack::DBSQL::ArchiveAdaptor") 
      unless($archive->isa("ReseqTrack::Archive"));
  my $sql = "delete from archive where archive_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $archive->dbID);
  $sth->execute;
  $sth->finish;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create an object for an undefined hashref") if(!$hashref);
  my $archive = ReseqTrack::Archive->new
      (
       -dbID => $hashref->{archive_id},
       -adaptor => $self,
       -md5 => $hashref->{md5},
       -relative_path => $hashref->{relative_path},
       -name => $hashref->{name},
       -file_id => $hashref->{file_id},
       -volume_name => $hashref->{volume_name},
       -created => $hashref->{created},
       -updated => $hashref->{updated},
       -priority => $hashref->{priority},
       -new_name => $hashref->{new_name},
       -new_relative_path => $hashref->{new_relative_path},
       -archive_action_id => $hashref->{archive_action_id},
       -archive_location_id => $hashref->{archive_location_id},
       -file_action_id => $hashref->{file_action_id},
       -file_exit_code => $hashref->{file_exit_code},
       -file_exit_reason => $hashref->{file_exit_reason},
      );
  return $archive;
}

sub delete_archive_lock{
  my ($self) = @_;
  #print "Deleting archive lock\n";
  my $ma = $self->db->get_MetaAdaptor;
  delete_lock_string("archive.lock", $ma);
}


1;
