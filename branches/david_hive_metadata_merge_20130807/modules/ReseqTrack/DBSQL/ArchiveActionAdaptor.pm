package ReseqTrack::DBSQL::ArchiveActionAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::ArchiveAction;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


sub columns{
  my ($self) = @_;
  return "archive_action.archive_action_id, archive_action.action";
}


sub table_name{
  my ($self) = @_;
  return "archive_action";
}

sub fetch_by_action{
  my ($self, $action) = @_;
  my $sql = "select ".$self->columns." from archive_action where action = ?";
  my $sth = $self->db->dbc->prepare($sql);
  $sth->bind_param(1, $action);
  $sth->execute();
  my $rowHashref = $sth->fetchrow_hashref;
  my $action_object = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;
  return $action_object;
}

sub store{
  my ($self, $archiveaction) = @_;
  my $sql = "insert into archive_action (archive_action_id, action) ".
      " values (?, ?)";
  my $sth = $self->db->dbc->prepare($sql);
  $sth->bind_param(1, $archiveaction->dbID);
  $sth->bind_param(2, $archiveaction->action);
  $sth->execute();
  $sth->finish;
  $archiveaction->adaptor($self);
  return $archiveaction;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create an object for an undefined hashref") if(!$hashref);
  my $action = ReseqTrack::ArchiveAction->new(
    -dbID => $hashref->{archive_action_id},
    -action => $hashref->{action},
    -adaptor => $self,
      );
  return $action;
}

1;
