package ReseqTrack::DBSQL::HostAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Host;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


#update

sub columns{
  return "host.host_id, host.name, host.remote, host.dropbox_dir";
}

sub table_name{
  return "host";
}


sub fetch_all_local{
  my ($self) = @_;
  my $sql = "select ".$self->columns." from host where remote = 0";
  my @hosts;
  my $sth = $self->prepare($sql);
  $sth->execute;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $host = $self->object_from_hashref($rowHashref);
    push(@hosts, $host);
  }
  $sth->finish;
  return \@hosts;
}


sub fetch_all_remote{
  my ($self) = @_;
  my $sql = "select ".$self->columns." from host where remote = 1";
  #print "sql is: $sql\n";
  my $sth = $self->prepare($sql);
  $sth->execute;
  my @hosts;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $host = $self->object_from_hashref($rowHashref) ;
#    print "host is " .$host->name . "\t host id is " . $host->dbID . "\n";
    push(@hosts, $host);
  }
  $sth->finish;
  return \@hosts;
}


sub fetch_by_name{
  my ($self, $name) = @_;
  my $sql = "select ".$self->columns." from host where name = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name);
  $sth->execute;
  my $rowHashref = $sth->fetchrow_hashref;
  my $host = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;
  return $host;
}


sub store{
  my ($self, $host) = @_;
  my $sql = "insert ignore into host ".
      "(name, remote) ".
      "values(?, ?) ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $host->name);
  $sth->bind_param(2, $host->remote);
  my $rows_inserted = $sth->execute();
  my $dbID = $sth->{'mysql_insertid'};
  if(!$dbID){
    my $new_host = $self->fetch_by_name($host->name);
    $host->dbID($new_host->dbID);
  }
  $sth->finish();
  $host->dbID($dbID) if($dbID);
  $host->adaptor($self);
}


sub update{
  my ($self, $host) = @_;
  throw("Can't update a host object without a dbID") unless($host->dbID);
  my $sql = "update host ".
      "set name = ?, ".
      "remote = ?".
     " where host_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $host->name);
  $sth->bind_param(2, $host->remote);
  $sth->bind_param(3, $host->dbID);
  $sth->execute();
  $sth->finish();
  return;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create an object from an undefined hashref") if(!$hashref);
  
  my $host = ReseqTrack::Host->new
      (
       -name => $hashref->{name},
       -adaptor => $self,
       -dbID => $hashref->{host_id},
       -remote => $hashref->{remote},
       -dropbox_dir => $hashref->{dropbox_dir},
      );
  return $host;
}



1;
