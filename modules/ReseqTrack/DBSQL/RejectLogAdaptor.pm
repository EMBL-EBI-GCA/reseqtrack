package ReseqTrack::DBSQL::RejectLogAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use lib '/homes/zheng/reseq-personal/zheng/lib/reseqtrack_hzb/modules/';

use ReseqTrack::RejectLog;
use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;
  my $self = $class->SUPER::new($db);
  return $self;
}

sub columns{
  return "reject_log.id, reject_log.file_id, reject_log.is_reject, ".
      "reject_log.reject_reason, reject_log.created";
}

sub table_name{
  return "reject_log";
}

sub fetch_by_file_id{
  my ($self, $file_id)  = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where file_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $file_id);
  $sth->execute;
  my @objects;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $object = $self->object_from_hashref($rowHashref) if($rowHashref);
    push(@objects, $object) if($object);
  }
  $sth->finish;
  return \@objects;
}

sub store{
  my ($self, $log_obj, $update) = @_;
  my $exists = $self->fetch_by_file_id($log_obj->file_id);  
  if($exists && @$exists == 1){
      warning("File already exists in reject_log table\n");
      if($update){
        $log_obj->dbID($exists->[0]->dbID);     
		return $self->update($log_obj);
      }   
  }
  elsif (@$exists > 1) {
      throw("More than 1 log objects exist for a file\n");
  }    
  
  throw("Can't store the file in reject_log table without an file id") unless($log_obj->file_id);
  
  my $sql = "insert into reject_log (file_id, is_reject, reject_reason, created) values(?, ?, ?, now())";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $log_obj->file_id);
  $sth->bind_param(2, $log_obj->is_reject);
  $sth->bind_param(3, $log_obj->reject_reason);
  $sth->execute;
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish;
  return $log_obj;
}

sub update{
  my ($self, $reject_log) = @_;

  my $sql = "update reject_log ".
      "set file_id = ? ".
      ", is_reject = ? ".
      ", reject_reason = ? ".
      ", created = now() ".
      "where id = ? ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $reject_log->file_id);
  $sth->bind_param(2, $reject_log->is_reject);
  $sth->bind_param(3, $reject_log->reject_reason);
  $sth->bind_param(4, $reject_log->dbID);  
  $sth->execute; 
  $sth->finish; 
  return $reject_log; 
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create a History object from an empty hashref") unless($hashref);
  my $reject_log = ReseqTrack::RejectLog->new(
    -dbID 			=>  $hashref->{id},
    -adaptor 		=> $self,
    -file_id 		=> $hashref->{file_id},
    -is_reject 		=> $hashref->{is_reject},
    -reject_reason 	=> $hashref->{reject_reason},
    -created 		=> $hashref->{created},
  );
  return $reject_log;
}

1;
