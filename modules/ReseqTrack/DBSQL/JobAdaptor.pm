package ReseqTrack::DBSQL::JobAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Job;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}


sub columns{
  return "job.job_id, job.submission_id, job.event_id, ".
      "job.stdout_file, job.stderr_file, job.input_string, job.exec_host, ".
      "job.retry_count, job_status.status, job_status.time ";
}

sub table_name{
  return "job, job_status";
}

sub where{
  return "job.job_id = job_status.job_id and job_status.is_current = 'y'";
}

sub internal_id_column{
  return "job.job_id";
}

sub fetch_by_input_string{
  my ($self, $dbID) = @_;
  my $sql = "select ".$self->columns.
      " from ".$self->table_name.
      " where ".$self->where.
      " and job.input_string = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $dbID);
  $sth->execute;
  my @jobs;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $job = $self->object_from_hashref($rowHashref);
    push(@jobs, $job);
  }
  $sth->finish;
  return \@jobs;
}

sub fetch_by_event_id{
  my ($self, $dbID) = @_;
  my $sql = "select ".$self->columns.
      " from ".$self->table_name.
      " where ".$self->where.
      " and job.event_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $dbID);
  $sth->execute;
  my @jobs;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $job = $self->object_from_hashref($rowHashref);
    push(@jobs, $job);
  }
  $sth->finish;
  return \@jobs;
}

sub fetch_by_event_name{
  my ($self, $event_name) = @_;
  my $aa = $self->db->get_EventAdaptor;
  my $event = $aa->fetch_by_event_name($event_name);
  throw("JobAdaptor:fetch_by_event_name, can't fetch jobs for ".$event_name.
        " as no associated event exists") if(!$event);
  return $self->fetch_by_event_id($event->dbID);
}

sub fetch_by_event{
  my ($self, $event) = @_;
  return $self->fetch_by_event_id($event->dbID);
}

sub fetch_by_current_status{
  my ($self, $status);
  my $sql = "select ".$self->columns.
      " from ".$self->table_name.
      " where ".$self->where.
      " and status = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $status);
  $sth->execute;
  my @jobs;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $job = $self->object_from_hashref($rowHashref);
    push(@jobs, $job);
  }
  $sth->finish;
  return \@jobs;
}


sub fetch_by_name_and_input_string{
  my ($self, $event_name, $input) = @_;
  my $sql = "select ".$self->columns.
      " from ".$self->table_name.", event".
      " where ".$self->where.     
      " and job.event_id = event.event_id ".
      " and event.name = ?".
      " and job.input_string = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $event_name);
  $sth->bind_param(2, $input);
  $sth->execute;
  my $rowHashref = $sth->fetchrow_hashref;
  my $job = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;
  return $job;
}

sub fetch_status{
  my ($self, $job) = @_;
  if(!$job->dbID){
    throw("Can't fetch status of ".$job." without a dbID");
  }
  my $sql = "select status from job_status where job_id = ? and is_current = 'y'";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $job->dbID);
  $sth->execute;
  my ($status) = $sth->fetchrow;
  $sth->finish;
  #throw("Failed to fetch status for ".$job->dbID) unless($status);
  $status = 'FINISHED' unless($status);
  $job->current_status($status) if($status);
  return $job;
}

sub store{
  my ($self, $job) = @_;
  my $job_insert_sql = "insert into job ".
      "(submission_id, event_id, input_string, stdout_file, stderr_file, ".
      "exec_host, retry_count)".
      "values(?, ?, ?, ?, ?, ?, ?)";
  throw("Can't store ".$job." without an event ".$job->event->name." dbID") 
      if(!$job->event->dbID);
  my $sth = $self->prepare($job_insert_sql);
  $sth->bind_param(1, $job->submission_id);
  $sth->bind_param(2, $job->event->dbID);
  $sth->bind_param(3, $job->input_string);
  $sth->bind_param(4, $job->stdout_file);
  $sth->bind_param(5, $job->stderr_file);
  $sth->bind_param(6, $job->host);
  $sth->bind_param(7, $job->retry_count);
  $sth->execute;
  my $dbID = $sth->{'mysql_insertid'};
  $sth->finish();
  $job->dbID($dbID),
  $job->adaptor($self);
  $job->current_status("CREATED") if(!$job->current_status);
  $self->set_status($job);
  return $job;
}

sub update{
  my ($self, $job) = @_;
  my $job_update_sql = "update job ".
      "set submission_id = ?, ".
      "event_id = ?, ".
      "input_string = ?, ".
      "stdout_file = ? ,".
      "stderr_file = ?,".
      "exec_host = ?, ".
      "retry_count = ? ".
      " where job_id = ? ";
  my $sth = $self->prepare($job_update_sql);
  $sth->bind_param(1, $job->submission_id);
  $sth->bind_param(2, $job->event->dbID);
  $sth->bind_param(3, $job->input_string);
  $sth->bind_param(4, $job->stdout_file);
  $sth->bind_param(5, $job->stderr_file);
  $sth->bind_param(6, $job->host);
  $sth->bind_param(7, $job->retry_count);
  $sth->bind_param(8, $job->dbID);
  $sth->execute;
  $sth->finish;
}

sub unset_submission_id{
  my ($self, $job) = @_;
  my $sql = "update job set submission_id = NULL";
  my $sth = $self->prepare($sql);
  $sth->execute;
  $sth->finish;
}
sub set_status{
  my ($self, $job) = @_;
  my $update_is_current = "update job_status ".
      "set is_current = 'n' ".
      "where is_current = 'y' ".
      "and job_id = ? ";
  my $sth = $self->prepare($update_is_current);
  $sth->bind_param(1, $job->dbID);
  $sth->execute;
  $sth->finish;
  my $job_status_insert = "insert into job_status ".
      "(job_id, status, time, is_current) ".
      "values(?, ?, now(), 'y')";
  $sth = $self->prepare($job_status_insert);
  $sth->bind_param(1, $job->dbID);
  $sth->bind_param(2, $job->current_status);
  $sth->execute;
  $sth->finish;
  return $job;
}


sub remove{
  my ($self, $job) = @_;
  print "Removing job ".$job->dbID."\n";
  my $sql = "delete from job where job_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $job->dbID);
  $sth->execute;
  $sth->finish;
  print "Removing job status ".$job->dbID."\n";
  my $status_sql = 'delete from job_status where job_id = ?';
  my $status_sth = $self->prepare($status_sql);
  $status_sth->bind_param(1, $job->dbID);
  $status_sth->execute;
  $status_sth->finish;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create an object from an undefined hashref") if(!$hashref);
  my $aa = $self->db->get_EventAdaptor;
  my $event = $aa->fetch_by_dbID($hashref->{event_id});
  warning("Couldn't find an event with event_id ".$hashref->{event_id})
      if(!$event);
  my $job = ReseqTrack::Job->new
      (
       -dbID => $hashref->{job_id},
       -adaptor => $self,
       -input_string => $hashref->{input_string},
       -submission_id => $hashref->{submission_id},
       -event => $event,
       -stdout => $hashref->{stdout_file},
       -stderr => $hashref->{stderr_file},
       -host => $hashref->{host},
       -current_status => $hashref->{status},
       -time => $hashref->{time},
       -retry_count => $hashref->{retry_count},
      );
  return $job;
}

;
