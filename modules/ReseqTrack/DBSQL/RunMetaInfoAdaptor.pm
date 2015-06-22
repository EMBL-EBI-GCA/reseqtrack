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
  return "run_meta_info.run_meta_info_id, run_meta_info.run_id, ".
      "run_meta_info.study_id, run_meta_info.study_name, run_meta_info.center_name, ".
      "run_meta_info.submission_id, run_meta_info.submission_date, ".
      "run_meta_info.sample_id, run_meta_info.sample_name, run_meta_info.population,".
      " run_meta_info.experiment_id, run_meta_info.instrument_platform, ".
      "run_meta_info.instrument_model, run_meta_info.library_name, ".
      "run_meta_info.run_name, run_meta_info.run_block_name, ".
      "run_meta_info.paired_length, run_meta_info.library_layout, ".
      "run_meta_info.status, run_meta_info.archive_base_count, ".
      "run_meta_info.archive_read_count, ".
	  "run_meta_info.library_strategy";
}

sub table_name{
  return "run_meta_info";
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
  return "insert into run_meta_info(run_id, study_id, study_name, center_name, ".
      "submission_id, submission_date, sample_id, sample_name, population, ".
      "experiment_id, instrument_platform, instrument_model, library_name, ".
      "run_name, run_block_name, paired_length, library_layout, status, ".
      "archive_base_count, archive_read_count, library_strategy) ".
      "values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
}

sub convert_date_sql{
   return "insert into run_meta_info(run_id, study_id, study_name, center_name, ".
      "submission_id, submission_date, sample_id, sample_name, population, ".
      "experiment_id, ".
      "instrument_platform, instrument_model, library_name, run_name, ".
      "run_block_name, paired_length, library_layout, status, archive_base_count, ".
      "archive_read_count, library_strategy) ".
      "values(?, ?, ?, ?, ?, STR_TO_DATE(?, '%Y-%M-%d %H:%i:%s'), ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
}
sub store{
  my ($self, $run_meta_info, $update, $convert_date) = @_;
  throw("Can't store without a run meta info object") unless($run_meta_info);
  my $existing = $self->fetch_by_run_id($run_meta_info->run_id);
  if($existing){
    $run_meta_info->dbID($existing->dbID);
    $run_meta_info->adaptor($self);
    if(are_run_meta_infos_identical($existing, $run_meta_info)){
      $self->store_history($run_meta_info);
      $self->store_statistics($run_meta_info, $update);
      return $run_meta_info;
    }else{
      if($update){
        return $self->update($run_meta_info);
      }else{
        warning("Can't update ".$run_meta_info." ".$run_meta_info->run_id.
                " already exists in the database and no update flag was passed");
        return undef;
      }
    }
  }
  my $sql = $self->store_sql;
  $sql = $self->convert_date_sql if($convert_date);
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $run_meta_info->run_id);
  $sth->bind_param(2, $run_meta_info->study_id);
  $sth->bind_param(3, $run_meta_info->study_name);
  $sth->bind_param(4, $run_meta_info->center_name);
  $sth->bind_param(5, $run_meta_info->submission_id);
  $sth->bind_param(6, $run_meta_info->submission_date);
  $sth->bind_param(7, $run_meta_info->sample_id);
  $sth->bind_param(8, $run_meta_info->sample_name);
  $sth->bind_param(9, $run_meta_info->population);
  $sth->bind_param(10, $run_meta_info->experiment_id);
  $sth->bind_param(11, $run_meta_info->instrument_platform);
  $sth->bind_param(12, $run_meta_info->instrument_model); 
  $sth->bind_param(13, $run_meta_info->library_name);
  $sth->bind_param(14, $run_meta_info->run_name);
  $sth->bind_param(15, $run_meta_info->run_block_name);
  $sth->bind_param(16, $run_meta_info->paired_length);
  $sth->bind_param(17, $run_meta_info->library_layout);
  $sth->bind_param(18, $run_meta_info->status);
  $sth->bind_param(19, $run_meta_info->archive_base_count);
  $sth->bind_param(20, $run_meta_info->archive_read_count);
  $sth->bind_param(21, $run_meta_info->library_strategy);	
  $sth->execute;
  my $dbID = $sth->{'mysql_insertid'};
  $run_meta_info->dbID($dbID);
  $run_meta_info->adaptor($self);
  $sth->finish;
  #storing any attached histories
  $self->store_history($run_meta_info);
  $self->store_statistics($run_meta_info, $update);
  #returning object with dbID and adaptor attached
  return $run_meta_info;
}

sub update{
  my ($self, $run_meta_info) = @_;
  throw("Can't update a ReseqTrack::RunMetaInfo object without a dbID") 
      unless($run_meta_info->dbID);
  my $existing = $self->fetch_by_dbID($run_meta_info->dbID);
  if(are_run_meta_infos_identical($run_meta_info, $existing)){
    return $run_meta_info;
  }
  if(!($run_meta_info->history) || @{$run_meta_info->history} == 0){
    throw("Can't change a run_meta_info object without  attaching a history object");
  }
  if(@{$run_meta_info->history} >= 1){
    my $no_dbID = 0;
    foreach my $history(@{$run_meta_info->history}){
      $no_dbID = 1 unless($history->dbID);
    }
    throw("All the history objects appear to already exist you need a new one")
        unless($no_dbID);
  }
  my $sql = "update run_meta_info ".
      "set study_id = ?, ".
      " study_name = ?, ".
      " center_name = ?, ".
      " submission_id = ?, ".
      " submission_date = ?, ".
      " sample_id = ?, ".
      " sample_name = ?, ".
      " population = ?, ".
      " experiment_id = ?, ".
      " instrument_platform = ?, ".
      " instrument_model = ?, ".
      " library_name = ?, ".
      " run_name = ?, ".
      " run_block_name = ?, ".
      " paired_length = ?, ".
      " library_layout = ?, ".
      " status = ?, ".
      " archive_base_count = ?, ".
      " archive_read_count = ?, ".
	  " library_strategy = ? ".
      "where run_meta_info_id = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $run_meta_info->study_id);
  $sth->bind_param(2, $run_meta_info->study_name);
  $sth->bind_param(3, $run_meta_info->center_name);
  $sth->bind_param(4, $run_meta_info->submission_id);
  $sth->bind_param(5, $run_meta_info->submission_date);
  $sth->bind_param(6, $run_meta_info->sample_id);
  $sth->bind_param(7, $run_meta_info->sample_name);
  $sth->bind_param(8, $run_meta_info->population);
  $sth->bind_param(9, $run_meta_info->experiment_id);
  $sth->bind_param(10, $run_meta_info->instrument_platform);
  $sth->bind_param(11, $run_meta_info->instrument_model); 
  $sth->bind_param(12, $run_meta_info->library_name);
  $sth->bind_param(13, $run_meta_info->run_name);
  $sth->bind_param(14, $run_meta_info->run_block_name);
  $sth->bind_param(15, $run_meta_info->paired_length);
  $sth->bind_param(16, $run_meta_info->library_layout);
  $sth->bind_param(17, $run_meta_info->status);
  $sth->bind_param(18, $run_meta_info->archive_base_count);
  $sth->bind_param(19, $run_meta_info->archive_read_count);
  $sth->bind_param(20, $run_meta_info->library_strategy);
  $sth->bind_param(21, $run_meta_info->dbID);
  $sth->execute;
  $sth->finish;
  my $hist_a = $self->db->get_HistoryAdaptor();
  $self->store_history($run_meta_info);
  $self->store_statistics($run_meta_info, 1);
  return $run_meta_info;
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


