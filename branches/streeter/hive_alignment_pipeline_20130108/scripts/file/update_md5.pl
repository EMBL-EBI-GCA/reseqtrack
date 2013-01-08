#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::IndexUtils qw(get_md5hash_from_sequence_index get_md5hash_from_alignment_index);
use ReseqTrack::DBSQL::DBAdaptor;

use ReseqTrack::File;
use ReseqTrack::Host;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $md5_file;
my $run;
my $help;
my $type;
my $skip_defined;
my $dir;

&GetOptions(
	    'dbhost=s'       => \$dbhost,
	    'dbname=s'       => \$dbname,
	    'dbuser=s'       => \$dbuser,
	    'dbpass=s'       => \$dbpass,
	    'dbport=s'       => \$dbport,
	    'md5_file=s' => \$md5_file,
	    'skip_defined!' => \$skip_defined,
	    'help!' => \$help,
	    'type=s' => \$type,
	    'dir=s' => \$dir,
	    'run!' => \$run,
    );

if($help){
  useage();
}


if(!$dbhost || !$dbname || !$dbuser){
  throw("Must provide database connection details with -dbhost -dbuser -dbpass ".
        "-dbport -dbname");
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );



my $md5_hash = get_md5hash($md5_file, 1);

my $fa = $db->get_FileAdaptor;

my $files;
if($type){
  $files = $fa->fetch_by_type($type);
}elsif($dir){
  my ($paths) = list_files_in_dir($dir, 1);
  foreach my $path(@$paths){
    my $file = $fa->fetch_by_name($path);
    unless($file){
      next;
    }
    push(@$files, $file) if($file);
  }
}else{
  $files = $fa->fetch_all;
}

my %file_hash;

foreach my $file(@$files){
  unless($file_hash{$file->name}){
    $file_hash{$file->filename} = $file;
  }else{
    print "Have duplicate for ".$file->name." in ".$file_hash{$file->filename}->name
        ."\n";
  }
}

foreach my $file(keys(%$md5_hash)){
  #my $name = basename($file);
  my $object = $file_hash{$file};
  if(!$object){
    #warning("Failed to find object for ".$file);
    next;
  }
  if($object->filename ne $file){
    warning($object->name." isn't the same as ".$file." skipping");
    next;
  }
  my $md5 = $md5_hash->{$file};
  if($skip_defined){
    next if($object->md5);
  }
  $object->md5($md5_hash->{$file});
  if($run){
    my $history = ReseqTrack::History->new(
      -other_id => $object->dbID,
      -table_name => 'file',
      -comment => "Updating md5 value",
        );
    $object->history($history);
    $fa->fast_update($object);
  }
}

sub useage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/file/update_md5.pl

=head1 SYNOPSIS

This script will use an md5 file to update the md5 value of the given files

=head1 OPTIONS

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the 
    standard port for mysql-g1kdcc.ebi.ac.uk
-md5_file, the file containing the md5s in the format md5sum filename
-type, To reduce database load the script fetches all the objects then updates the
ones it has md5s for. If the type is specified, then only objects with this type
are fetched so the set to be updated must all be of this typ
-run, Binary flag to indicate sql should be executed
-help, Binary flag to indicate if the help should be printed out

=head1 Examples

This would update all the files in the md5 file regardless of the type

perl ReseqTrack/scripts/file/update_md5.pl -dbhost myhost -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname my_database -md5_file files.md5

will only update those objects which have type FASTQ

perl ReseqTrack/scripts/file/update_md5.pl -dbhost myhost -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname my_database -md5_file files.md5 -type FASTQ

=cut

