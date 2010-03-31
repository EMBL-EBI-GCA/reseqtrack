#!/usr/local/bin/perl -w

use strict;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::ArchiveUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::GeneralUtils;
use File::Basename;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my @files;
my $file_list;
my $dir;
my $run = 0;
my $descend = 1;
my $from_db;
my $type;
my $path_like;
my $action_string;
my $action_location_name;
my $sleep = 240;
my $skip_cleanup = 0;
my $verbose = 0;
my $create_changelog_detailed = 1;
my $changelog_path = '/nfs/1000g-work/G1K/archive_staging/ftp/changelog_details/';
my $changelog_name;
my $update_changelog = 1;
my $original_changelog = '/nfs/1000g-archive/vol1/ftp/CHANGELOG';
my $new_changelog = '/nfs/1000g-work/G1K/archive_staging/ftp/CHANGELOG';
my $priority = 50;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'file=s@' => \@files,
  'file_list=s' => \$file_list,
  'dir=s' => \$dir,
  'descend!' => \$descend,
  'run!' => \$run,
  'from_db!' => \$from_db,
  'type=s' => \$type,
  'path_like=s' => \$path_like,
  'action=s' => \$action_string,
  'sleep=s' => \$sleep,
  'skip_cleanup!' => \$skip_cleanup,
  'verbose!' => \$verbose,
  'create_changelog!' => \$create_changelog_detailed,
  'changelog_path=s' => \$changelog_path,
  'changelog_name=s' => \$changelog_name,
  'update_changelog!' => \$update_changelog,
  'original_changelog=s' => \$original_changelog,
  'new_changelog=s' => \$new_changelog,
  'priority=s' => \$priority,
    );

my $date = current_time;
my @values = split /\s+/, $date;
my $day = $values[0];
my $timestamp = $day;
$timestamp =~ s/\-//g;

#connect to db
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

my $archive_action_adaptor = $db->get_ArchiveActionAdaptor;
my $archive_location_adaptor = $db->get_ArchiveLocationAdaptor;

my $archive_actions = $archive_action_adaptor->fetch_all;
my %action_hash;

foreach my $action(@$archive_actions){
  $action_hash{$action->action} = $action;
}

my $other_location_name;

if($action_string){
  if($action_string eq "archive" || $action_string eq 'replace'){
    $action_location_name = "staging";
    $other_location_name = "archive";
  }elsif($action_string eq "dearchive"){
    $action_location_name = "archive";
    $other_location_name = "staging";
  }else{
    throw("Don't know what to do with action ".$action_string.
          " should be archive or dearchive");
  }
}




my $archive_location = $archive_location_adaptor->fetch_by_archive_location_name($action_location_name);
my $other_location = $archive_location_adaptor->fetch_by_archive_location_name($other_location_name);
my $location_root = $archive_location->location;
my $new_root = $other_location->location;
if($new_root =~ /\/nfs\/1000g-archive/){
  $new_root .= "/vol1";
}
my $aa = $db->get_ArchiveAdaptor;
my $archives = $aa->fetch_all;
cleanup_archive($archives, $db, 0);

my $fa = $db->get_FileAdaptor;  
#generate file list
my %file_objects;
if($file_list){
  my $list = get_lines_from_file($file_list);
  push(@files, @$list);
}elsif($dir){
  my ($list, $hash) = list_files_in_dir($dir, 1);
  if($descend){
    push(@files, @$list);
  }else{
    my $dir_list = $hash->{$dir};
    push(@files, @$dir_list);
  }
}elsif($from_db){
  my $list;
  if($path_like && !$type){
    $list = $fa->fetch_all_like_path($path_like);
  }elsif(!$path_like && $type){
    $list = $fa->fetch_by_type($type);
  }elsif($path_like && $type){
    my $temps = $fa->fetch_all_like_path($path_like);
    foreach my $temp(@$temps){
      next unless($temp->type eq $type);
      push(@$list, $temp);
    }
  }else{
   print "When archiving files from the database you must specify a type and/or ".
       "a root path for the files to be from\n";
  }
  print "Have ".@$list." files to check\n" if($verbose);
  foreach my $file(@$list){
    my $other_root = $other_location->location;
    next if($file->path =~ /$other_root/);
    push(@files, $file->full_path);
    $file_objects{$file->full_path} = $file;
  }
}

throw("Need more than zero froms in file array from either the -file, ".
      "-file_list or -dir") unless(@files >= 1);

#sanity check files
my @files_to_archive;
my %which_action_hash;
my %changelog_hash;
foreach my $file(@files){
  warning("Can't archive a file ".$file." which doesn't exist") unless(-e $file);
  unless($file =~ /$location_root/){
    throw("Can't ".$action_string." ".$file." it isn't in the ".$location_root);
  }
  my $new_file = $file;
  $new_file =~ s/$location_root/$new_root/;
  if($archive_location->location_name eq 'staging'){
    if(-e $new_file){
      #print STDERR "Running replace on ".$file."\n" 
      #    unless($action_string eq 'replace');
      $which_action_hash{$file} = $action_hash{'replace'};
    }else{
      $which_action_hash{$file} = $action_hash{'archive'};
    }
  }else{
    if(-e $new_file){
      warning("Can't dearchive ".$file." as ".$new_file." exists");
      next;
    }else{
      $which_action_hash{$file} = $action_hash{'dearchive'};
    }
  }
  push(@files_to_archive, $file);
}


foreach my $file_path(@files_to_archive){
  next unless(-e $file_path);
  my $action = $which_action_hash{$file_path};
  $file_path =~ s/\/$//;
  my $file;
  unless($changelog_hash{$action->action}){
    $changelog_hash{$action->action} = [];
  }
  unless($file_objects{$file_path}){
    $file = $fa->fetch_by_name($file_path);
    if(!$file){
      throw("Failed to fetch file from ".$file_path);
    }
  }else{
    $file = $file_objects{$file_path};
  }
  push(@{$changelog_hash{$action->action}}, $file_path);
  $changelog_name = lc($file->type) unless($changelog_name);
  my $archive = create_archive_from_objects($file, $action, $archive_location);
  $archive->priority($priority);
  $aa->store($archive) if($run);
}


unless($skip_cleanup){
  my $find_all = 1;
 UPDATE:while($find_all){
   my $archives = $aa->fetch_all;
   $find_all = 0 if(!$archives || @$archives == 0);
   cleanup_archive($archives, $db, $verbose);
   if(!$find_all){
     last UPDATE;
   }
   print "Sleep $sleep\n";
   sleep($sleep);
 }
}

my @changelog_files;
if($create_changelog_detailed){
   my $root_file_path = $changelog_path."/changelog_details_".$timestamp;
  foreach my $action(keys(%changelog_hash)){
     my $type = get_changelog_name($action);
    my $full_path = $root_file_path."_".$type."_".$changelog_name;
    $full_path =~ s/\/\//\//g;
     throw("Can't write ".$full_path." as it already exists") if(-e $full_path);
    if($run){
      open(FH, ">".$full_path) or throw("Failed to open ".$full_path." $!");
      foreach my $file(@{$changelog_hash{$action}}){
        $file =~ s/\/nfs\/1000g\-//;
        $file =~ s/archive\///;
        $file =~ s/work\/G1K\///;
        $file =~ s/archive\_staging\///;
        $file =~ s/vol1\///;
        $file =~ s/ftp\///;
        print FH $file."\n";
      }
      close(FH);
    }
    push(@changelog_files, $full_path);
  }
}

if($update_changelog){
   if(@changelog_files >= 1 && $run){
    open(FH, $original_changelog) or throw("Failed to open ".$original_changelog." $!");
    my @data = <FH>;
    close(FH);
    open(FH, ">",$new_changelog) or throw("Failed to open ".$new_changelog." $! ");
    print FH "\n".$day."\n\n";
    print FH "Modification to ".$changelog_name."\n";
    print FH "Details can be found in the changelog details files\n\n";
    my $trim_root = $changelog_path;
    $trim_root =~ s/changelog_details\///;
    foreach my $file(@changelog_files){
      $file =~ s/$trim_root//;
      print FH $file."\n";
    }
    print FH "\n";
    foreach my $line(@data){
      print FH $line;
    }
    close(FH);
  }
}
#print "Deleting archive lock\n";
$aa->delete_archive_lock;
#print "Deleted archive lock\n";
sub get_changelog_name{
  my ($action) = @_;
  if($action eq 'archive'){
    return "new";
  }elsif($action eq 'replace'){
    return "replacement";
  }elsif($action eq 'dearchive'){
    return "withdrawn";
  }
}



