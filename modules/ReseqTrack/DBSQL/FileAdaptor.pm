package ReseqTrack::DBSQL::FileAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::File;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileUtils qw(are_files_identical);
use File::Basename;

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub columns{
  return "file.file_id, file.name, file.md5, file.type, file.size, file.host_id, ".
      "file.withdrawn, file.created, file.updated"
}

sub table_name{
  return "file";
}

sub fetch_by_name{
  my ($self, $name) = @_;
  my $sql = "select ".$self->columns." from file ".
      "where name = ? ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name);
  $sth->execute;
  my @files;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $file = $self->object_from_hashref($rowHashref);
    push(@files, $file);
  }
  $sth->finish;
  throw($name." has returned multiple objects ".@files." not sure what to do") 
      if(@files && @files >= 2);
  my $file = $files[0];
  return $file;
}


sub fetch_by_name_and_md5{
  my ($self, $name, $md5) = @_;
  my $sql = "select ".$self->columns." from file ".
      "where name = ?  and md5 = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name);
  $sth->bind_param(2, $md5);
  $sth->execute;
  my $rowHashref = $sth->fetchrow_hashref;
  my $file = $self->object_from_hashref($rowHashref) if($rowHashref);
  $sth->finish;
  return $file;
}

sub fetch_by_type{
  my ($self, $type) = @_;
  my $sql = "select ".$self->columns." from file ".
      "where type = ? ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $type);
  $sth->execute;
  my @files;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $file = $self->object_from_hashref($rowHashref);
    push(@files, $file);
  }
  $sth->finish;
  return \@files;
}

sub fetch_by_filename{
  my ($self, $name) = @_;
  my $sql = "select ".$self->columns." from file ".
      "where name like ? ";
  my $name_like = '%/'.$name;
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name_like);
  $sth->execute;
  my @files;
  while(my $rowHashref = $sth->fetchrow_hashref){
    #print "Have rowHashref ".$rowHashref."\n";
    my $file = $self->object_from_hashref($rowHashref);
    #print "Getting back ".$file."\n";
    push(@files, $file);
  }
  $sth->finish;
  #print "Returning @files ".@files."\n";
  #print $files[0]."\n";
  return \@files;
}

sub fetch_by_dirname{
  my ($self, $name) = @_;
  my $sql = "select ".$self->columns." from file ".
      "where name like ? ";
  my $name_like = $name."%";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name_like);
  $sth->execute;
  my @files;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $file = $self->object_from_hashref($rowHashref);
    push(@files, $file);
  }
  $sth->finish;
  return \@files;
}

sub fetch_all_like_name{
  my ($self, $name) = @_;
  my $sql = "select ".$self->columns." from file ".
      "where name like ? ";
  my $name_like = "%".$name."%";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name_like);
  $sth->execute;
  my @files;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $file = $self->object_from_hashref($rowHashref);
    push(@files, $file);
  }
  $sth->finish;
  return \@files;
}


sub fetch_all_like_path{
   my ($self, $name) = @_;
  my $sql = "select ".$self->columns." from file ".
      "where name like ? ";
  my $name_like = $name."%";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $name_like);
  $sth->execute;
  my @files;
  while(my $rowHashref = $sth->fetchrow_hashref){
    my $file = $self->object_from_hashref($rowHashref);
    push(@files, $file);
  }
  $sth->finish;
  return \@files;
}

sub fetch_by_host {
	my ($self, $host_id) = @_;
	my $sql = "select " .$self->columns. " from file ".
		"where host_id = ? ";
	my $sth = $self->prepare($sql);
	$sth->bind_param(1, $host_id);
	$sth->execute; 
	my @files;
	while(my $rowHashref = $sth->fetchrow_hashref){
    	my $file = $self->object_from_hashref($rowHashref);
    	push(@files, $file);
	}
	$sth->finish;
	return \@files;
}	


sub store{
  my ($self, $file, $update, $dont) = @_;
  #check if file object already exists
  #does the exact path exist
  #does something with the same name exist
  unless($dont){
    my $existing = $self->fetch_by_name($file->name);
    if($existing){
      $file->dbID($existing->dbID);
      $file->adaptor($self);
      $file->created($existing->created);
      if(are_files_identical($existing, $file)){
        $self->store_history($file);
        $self->store_statistics($file, $update);
        return $file;
      }
      if($update){
        return $self->update($file);
      }else{
        warning($file." ".$file->name." wasn't stored as the name already exists ".
                " in the database and you are disallowing update");
        return undef;
      }
    }else{
      my $possible_existing = $self->fetch_by_filename($file->filename);
      if($possible_existing){
        if(@$possible_existing == 1){
          my $existing = $possible_existing->[0];
          $file->dbID($existing->dbID);
          $file->adaptor($self);
          $file->created($existing->created);
          if(are_files_identical($existing, $file)){
            $self->store_history($file);
            $self->store_statistics($file, $update);
            return $file;
          }
          if($update){
            return $self->update($file);
          }else{
            warning($file." ".$file->name." wasn't stored as the ".
                    "filename already exists  in the database and you are ".
                    "disallowing update");
            return undef;
          }
        }elsif(@$possible_existing >= 2){
          my $for_update;
          foreach my $existing(@$possible_existing){
            if($existing->type eq $file->type){
              if($for_update){
                warning("Can't update ".$file->filename." there are multiple files ".
                        "which share its name and type");
              }
              $for_update = $existing;
            }
          } 
          $file->dbID($for_update->dbID) if($for_update);
          return $self->update($file) if($for_update);
          warning("There are multiple files in the database with the filename ".
                  $file->filename." and non with a matching type ".$file->type.
                  " Don't know what to do now");
          return undef;
        }
      }
    }
  }
  #removing any the final / in paths
  if($file->name =~ /\/$/){
    my $tmp = $file->name;
    $tmp =~ s/\/$//;
    $file->name($tmp);
  }
  #removing any double slashes
  if($file->name =~ /\/\//){
    my $tmp = $file->name;
    $tmp =~ s/\/\//\//g;
    $file->name($tmp);
  }
  #checking to see if the hose needs to be stored
  if(!$file->host->dbID){
    my $ha = $self->db->get_HostAdaptor;
    $ha->store($file->host);
  }
  throw("Can't store ".$file." as doesn't have a host dbID") unless($file->host->dbID);
  #sql to store the object
  my $sql = "insert into file (name, md5, type, size, host_id, withdrawn, ".
      "created, updated) values(?, ?, ?, ?, ?, ?, now(), now())";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $file->name);
  $sth->bind_param(2, $file->md5);
  $sth->bind_param(3, $file->type);
  $sth->bind_param(4, $file->size);
  $sth->bind_param(5, $file->host->dbID);
  $sth->bind_param(6, $file->withdrawn);
  $sth->execute;
  my $dbID = $sth->{'mysql_insertid'};
  $file->dbID($dbID);
  $file->adaptor($self);
  $sth->finish;
  #storing any attached histories
  $self->store_history($file);
  $self->store_statistics($file);
  #returning object with dbID and adaptor attached
  return $file;
}

sub update_withdrawn{
  my ($self, $file) = @_;
  $self->update($file, 1);
}

sub update{
  my ($self, $file, $allow_withdraw, $change_name) = @_;
  my $original;
  if(!$file->dbID){
    $original = $self->fetch_by_name_and_md5($file->name, $file->md5);
    if($original){
      $file->dbID($original->dbID);
      throw($original." ".$original->dbID." has no created timestamp")
          unless($original->created);
      $file->created($original->created);
    }else{
      $original = $self->fetch_by_name($file->name);
      if($original){
        $file->dbID($original->dbID);
         throw($original." ".$original->dbID." has no created timestamp")
             unless($original->created);
        $file->created($original->created);
      }else{
        my $others = $self->fetch_by_filename($file->filename);
        if($others && @$others == 1){
          $file->dbID($others->[0]->dbID);
           throw($others->[0]." ".$others->[0]->dbID." has no created timestamp")
               unless($others->[0]->created);
          $file->created($others->[0]->created);
        }else{
          throw($file->full_path." ".$file->md5." doesn't yet existing in the ".
                "database you can't update it");
        }
      }
    }
  }
  $original = $self->fetch_by_dbID($file->dbID);
  unless($file->created){
    $file->created($original->created);
  }
  if($original->filename ne $file->filename){
    throw("Can't change the name of a file from ".$original->name." to ".$file->name.
          " using update unless change_name is ".
          "specifed") unless($change_name);
  }
  if($original->withdrawn != $file->withdrawn){
    throw("Can't alter withdrawal status of ".$file->full_path.
          " unless allow_withdraw is set") unless($allow_withdraw);
  }
  unless($original->filename ne $file->filename){
    if(are_files_identical($original, $file)){
      throw("Not updating ".$file->dbID." ".$file->name." it is exactly the same as ".
            "the original");
    }
  }
  if(!($file->history) || @{$file->history} == 0){
    throw("Can't change a file object without  attaching a history object");
  }
  if(@{$file->history} >= 1){
    my $no_dbID = 0;
    foreach my $history(@{$file->history}){
      $no_dbID = 1 unless($history->dbID);
    }
    throw("All the history objects appear to already exist you need a new one")
        unless($no_dbID);
  }
  if(!$file->host->dbID){
    my $ha = $self->db->get_HostAdaptor;
    $ha->store($file->host);
  }
  #removing any the final / in paths
  if($file->name =~ /\/$/){
    my $tmp = $file->name;
    $tmp =~ s/\/$//;
    $file->name($tmp);
  }
  #removing any double slashes
  if($file->name =~ /\/\//){
    my $tmp = $file->name;
    $tmp =~ s/\/\//\//g;
    $file->name($tmp);
  }
  unless($file->created){
    throw("Can't update ".$file->name." ".$file->dbID." it has no created timestamp");
  }
  my $sql = "update file ".
      "set md5 = ?, ".
      "name = ?, ".
      "type = ?, ".
      "size = ?, ".
      "host_id = ?, ".
      "withdrawn = ?, ".
      "created = ?, ".
      "updated = now() ".
      "where file_id = ? ";
  my $sth = $self->prepare($sql);
  
  $sth->bind_param(1, $file->md5);
  $sth->bind_param(2, $file->name);
  $sth->bind_param(3, $file->type);
  $sth->bind_param(4, $file->size);
  $sth->bind_param(5, $file->host->dbID);
  $sth->bind_param(6, $file->withdrawn);
  $sth->bind_param(7, $file->created);
  $sth->bind_param(8, $file->dbID);
  $sth->execute;
  $sth->finish;
  #storing any attached histories
  $self->store_history($file);
  $self->store_statistics($file, 1);
  return $file;
}

sub fast_update{
  my ($self, $file) = @_;
  if(!$file->dbID){
    throw("Can't call ReseqTrack::DBSQL::FileAdaptor->fast_update without an dbID ".
	  "already associated with your ".$file);
  }
  throw("Can't update a file if created isn't already already defined")
    unless($file->created);
  #removing any double slashes
  if($file->name =~ /\/\//){
    my $tmp = $file->name;
    $tmp =~ s/\/\//\//g;
    $file->name($tmp);
  }
  my $sql = "update file ".
      "set md5 = ?, ".
      "name = ?, ".
      "type = ?, ".
      "size = ?, ".
      "host_id = ?, ".
      "withdrawn = ?, ".
      "created = ?, ".
      "updated = now() ".
      "where file_id = ? ";
  my $sth = $self->prepare($sql);
  
  $sth->bind_param(1, $file->md5);
  $sth->bind_param(2, $file->name);
  $sth->bind_param(3, $file->type);
  $sth->bind_param(4, $file->size);
  $sth->bind_param(5, $file->host->dbID);
  $sth->bind_param(6, $file->withdrawn);
  $sth->bind_param(7, $file->created);
  $sth->bind_param(8, $file->dbID);
  $sth->execute;
  $sth->finish;
  #storing any attached histories
  $self->store_history($file);
  $self->store_statistics($file, 1);
  return $file;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create a ReseqTrack::File from an undefined hashref") if(!$hashref);
  my $file = ReseqTrack::File->new(
				   -host_id => $hashref->{host_id},
				   -dbID => $hashref->{file_id},
				   -adaptor => $self,
				   -name => $hashref->{name},
				   -md5 => $hashref->{md5},
				   -type => $hashref->{type},
				   -size => $hashref->{size},
				   -path => $hashref->{path},
				   -withdrawn => $hashref->{withdrawn},
				   -created => $hashref->{created},
				   -updated => $hashref->{updated}
				  );
  return $file;
}

sub change_file_type {
  my ($self, $file, $new_type) = @_;
  throw ("new type not specified") if (! $new_type);

  my $fha = $self->db->get_FileHistoryAdaptor;
  my $fh = OneKGenomes::FileTracking::FileHistory->new
      (
       -first_file_id => $file->dbID,
       -second_file_id=> undef,
       -comment=>'type changed',
       -adaptor=>$fha,
      );
  my $change_is_current_sql = "update file set type = ? where unique_ident_id = ?";
  my $change_sth = $self->prepare($change_is_current_sql);
  $change_sth->bind_param(1, $new_type);
  $change_sth->bind_param(2, $file->unique_ident->dbID);
  $change_sth->execute;
  $change_sth->finish;
  $fha->store($fh);
}


1;
