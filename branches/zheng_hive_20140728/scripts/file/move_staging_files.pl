#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::Exception;
use File::Basename;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $move_list;
my $run = 0;
my $from;
my $to;
my $clobber = 0;

&GetOptions( 
	    'dbhost=s'      => \$dbhost,
	    'dbname=s'      => \$dbname,
	    'dbuser=s'      => \$dbuser,
	    'dbpass=s'      => \$dbpass,
	    'dbport=s'      => \$dbport,
	    'move_list=s' => \$move_list,
	    'from=s' => \$from,
	    'to=s' => \$to,
	    'run!' => \$run,
	    'clobber!' => \$clobber,
    );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

my $hash;
if($move_list){
  $hash = parse_movelist($move_list);
}elsif($from && $to){ 
  $hash->{$from} = $to;
}else{
  throw("Need to give move_archive_files.pl either a -move_list or a -from and ".
	"a -to value");
}
my $fa = $db->get_FileAdaptor;

my @file_objects;
foreach my $path(keys(%$hash)){
  print $path."\n";
  my $new_path = $hash->{$path};
  unless($clobber){
    next if(-e $new_path);
  }
  my $file = $fa->fetch_by_name($path);
  unless($file){
    throw("Failed to find a file object for ".$path);
  }
  my $new_dir = dirname($new_path);
  my $new_name = basename($new_path);
#  unless($new_name eq $file->filename){
#    throw("Can't rename files using this script".$path." to ".$new_path.
#	  " is not possible");
#  }
  move_file_in_db_and_dir([$file], $new_dir, $file->type, $db);
  unless(-e $new_path){
    throw("Failed to move ".$path." to ".$new_path);
  }
}
