package ReseqTrack::DBSQL::RunMetaInfoAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::RunMetaInfo;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunMetaInfoUtils qw (are_run_meta_infos_identical);

use File::Basename;

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}



sub columns{
  my ($self) = @_;
  my $table_name = $self->table_name;
    return join(', ', map {"$table_name.$_"}
      qw(run_meta_info_id run_id study_id study_name center_name submission_id
      submission_date sample_id sample_name population experiment_id
      instrument_platform instrument_model library_name run_name run_block_name
      paired_length library_layout status archive_base_count archive_read_count
      library_strategy)
      );
}

sub table_name{
  return "run_meta_info_vw";
}

sub fetch_by_name{
  my ($self, $name) = @_;
  return $self->fetch_by_run_id($name);
}

sub fetch_by_run_id{
  my ($self, $run_id) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where run_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $run_id);
  $sth->execute;
  my $rowhashref = $sth->fetchrow_hashref;
  my $run_meta_info = $self->object_from_hashref($rowhashref) if($rowhashref);
  $sth->finish;
  return $run_meta_info;
 
}

sub fetch_incomplete_by_event{
  my ($self, $event) = @_;
  throw($event->name . " table name is not " . $self->table_name) if ($event->table_name ne $self->table_name);
  my $table_name = $self->table_name;
  my $sql = "select ".$self->columns." from $table_name ".
      "left outer join (select event_complete.other_id ".
      "from event_complete where event_complete.event_id = ?) as e ".
      "on $table_name.run_meta_info_id = e.other_id ".
      "where e.other_id is null";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $event->dbID);
  $sth->execute;
  my @run_meta_infos;
  while(my $hashref = $sth->fetchrow_hashref){
    my $run_meta_info = $self->object_from_hashref($hashref);
    push(@run_meta_infos, $run_meta_info);
  }
  
  $sth->finish;
  return \@run_meta_infos;
}

sub fetch_by_sample_id{
  my ($self, $sample_id) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where sample_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $sample_id);
  $sth->execute;
  my @run_meta_infos;
  while(my $rowhashref = $sth->fetchrow_hashref){
    my $run_meta_info = $self->object_from_hashref($rowhashref) if($rowhashref);
    push(@run_meta_infos, $run_meta_info) if($run_meta_info);
  }
  $sth->finish;
  return \@run_meta_infos;
}

sub fetch_by_sample_name{
  my ($self, $sample_name) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where sample_name = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $sample_name);
  $sth->execute;
  my @run_meta_infos;
  while(my $rowhashref = $sth->fetchrow_hashref){
    my $run_meta_info = $self->object_from_hashref($rowhashref) if($rowhashref);
    push(@run_meta_infos, $run_meta_info) if($run_meta_info);
  }
  $sth->finish;
  return \@run_meta_infos;
}

sub fetch_by_study_id{
  my ($self, $study_id) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where study_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $study_id);
  $sth->execute;
  my @run_meta_infos;
  while(my $rowhashref = $sth->fetchrow_hashref){
    my $run_meta_info = $self->object_from_hashref($rowhashref) if($rowhashref);
    push(@run_meta_infos, $run_meta_info) if($run_meta_info);
  }
  $sth->finish;
  return \@run_meta_infos;
}

sub fetch_by_submission_id{
  my ($self, $submission_id) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.
      " where submission_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $submission_id);
  $sth->execute;
  my @run_meta_infos;
  while(my $rowhashref = $sth->fetchrow_hashref){
    my $run_meta_info = $self->object_from_hashref($rowhashref) if($rowhashref);
    push(@run_meta_infos, $run_meta_info) if($run_meta_info);
  }
  $sth->finish;
  return \@run_meta_infos;
}

sub store_sql{
  throw("Cannot store - run_meta_info has been superseded");
}

sub convert_date_sql{
  throw("Cannot store - run_meta_info has been superseded");
}
sub store{
  throw("Cannot store - run_meta_info has been superseded");
}

sub update{
  throw("Cannot update - run_meta_info has been superseded");
}



sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create an object from an empty hashref") 
      unless($hashref && keys(%$hashref) >= 1);
  my $run_meta_info = ReseqTrack::RunMetaInfo->new(
    -dbID => $hashref->{run_meta_info_id},
    -adaptor => $self,
    -run_id => $hashref->{run_id},
    -study_id => $hashref->{study_id},
    -study_name => $hashref->{study_name},
    -center_name => $hashref->{center_name},
    -submission_id => $hashref->{submission_id},
    -submission_date => $hashref->{submission_date},
    -sample_id => $hashref->{sample_id},
    -sample_name => $hashref->{sample_name},
    -population => $hashref->{population},
    -experiment_id => $hashref->{experiment_id},
    -instrument_platform => $hashref->{instrument_platform},
    -instrument_model => $hashref->{instrument_model},
    -library_name => $hashref->{library_name},
    -run_name => $hashref->{run_name},
    -run_block_name => $hashref->{run_block_name},
    -paired_length => $hashref->{paired_length},
    -library_layout => $hashref->{library_layout},
    -status => $hashref->{status},
    -archive_base_count => $hashref->{archive_base_count},
    -archive_read_count => $hashref->{archive_read_count},
	-library_strategy => $hashref->{library_strategy},
      );
  return $run_meta_info;
}


