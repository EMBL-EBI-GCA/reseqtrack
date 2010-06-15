package ReseqTrack::Tools::ArchiveUtils;

use strict;
use warnings;
use Exporter;
use File::Basename;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Archive;
use ReseqTrack::File;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(create_archive_from_objects create_file_from_archive cleanup_archive);


sub create_archive_from_objects{
  my ($file, $archive_action, $archive_location, $new_path) = @_;
  throw("Can't archive without a file object") unless($file);
  throw("Shouldn't be archiving a file which doesn't exist ".$file->full_path)
      unless(-e $file->full_path);
  my $archive = ReseqTrack::Archive->new
      (
       -file => $file,
       -archive_action => $archive_action,
       -archive_location => $archive_location,
      );
  if($new_path){
    my ($new_name, $new_dir) = fileparse($new_path);
    $archive->new_name($new_name);
    if($archive->name ne $archive->new_name){
      warning("Changing the name of ".$file->full_path." to ".$archive->new_name);
    }
    my $root = $archive_location->location;
    my $relative_path = $new_dir;
    $relative_path =~ s/$root//;
    $relative_path =~ s/\/vol\d+//;
    $archive->new_relative_path($relative_path);
  }
  if(!$archive->relative_path){
    throw("Can't archive ".$archive->full_path." don't have a relative path");
  }
  return $archive;
}

sub create_file_from_archive{
  my ($archive, $db, $host, $location) = @_;
  #print "Creating file from ".$archive."\n";
  my $full_path = $archive->full_path;
  my $filename = basename($full_path);
  my $dirname = dirname($full_path);
  unless(-e $full_path){
    throw("Can't create a file for ".$full_path." it doesn't exist");
    #return undef;
  }
  throw("Can't create file from archive without a host object") 
      unless($host && $host->isa("ReseqTrack::Host"));
  throw("Can't create a file from archive without a dbadaptor")
      unless($db && $db->isa("ReseqTrack::DBSQL::DBAdaptor"));
  #try and find ident
  #print "Trying to fetch ".$archive->file_id." from ".$db->dbc->dbname."\n";
  my $current_file = $db->get_FileAdaptor->fetch_by_dbID($archive->file_id);
  
  if(!$archive->md5){
    print "Problem we seem to have no md5\n";
  }
  
  #print "Creating file object\n";
  my $file_object = ReseqTrack::File->new
      (
       -name => $full_path,
       -type => $current_file->type,
       -md5 => $archive->md5,
       -host => $host,
      );
  #print "Have file object with ".$file_object->full_path."\n";
  return $file_object;
}


sub cleanup_archive{
  my ($archives, $db, $verbose) = @_;
  my $fa = $db->get_FileAdaptor;
  my $aa = $db->get_ArchiveAdaptor;
 ARCHIVE:foreach my $archive(@$archives){
    my $old_file = $fa->fetch_by_dbID($archive->file_id);
    if(!$old_file){
      throw("Can't run ".$archive->name." ".$archive->dbID." doesn't have an ".
            "associated file");
    }
    my $old_path = $old_file->full_path;
    my $new_path = $archive->full_path;
    $old_path =~ s/\/\//\//g;
    $new_path =~ s/\/\//\//g;
    if($verbose){
      print "Old path ".$old_path."\n";
      print "New path ".$new_path."\n";
    }
    if($old_path eq $new_path){
      #print "Skipping ".$archive->dbID." old path ".$old_path." eq ".$new_path."\n";
      next ARCHIVE;
    }
    my $new_file;
    eval{
      #print "Trying to create new file from ".$archive->dbID."\n";
      $new_file = create_file_from_archive($archive, $db, $old_file->host);
      unless($new_file){
        throw("Seem to of failed to create a new file object from ".$archive->dbID);
      }
    };
    if($@){
      warning("Problem with archiving ".$archive->file_id."\n $@");
      next ARCHIVE;
    }
    if($new_file && -e $new_file->full_path){
      my $history;
      if($new_file->filename eq $old_file->filename){
        $history = create_history($new_file, $old_file, $old_path." mv to ".$new_path);
      }else{
        $history = create_history(undef, $old_file, "name changed from ".$old_file->name." to ".$new_file->name);
      }
      $new_file->history($history);
      $new_file->dbID($old_file->dbID);
      $new_file->created($old_file->created);
      my $return = $fa->fast_update($new_file);
      if($return){
        $aa->remove($archive);
      }else{
        throw("Failed to store ".$new_file->name);
      }
    }else{
      print STDERR $new_file->full_path."\n";
      print STDERR "***PROBLEMS WITH ABOVE PATH, it doesn't exist****\n";
    }
  }
}
