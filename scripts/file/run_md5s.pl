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
my $filename;
my $run;
my $help;
my $type;
my $check_all = 0;
&GetOptions(
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'file=s' => \$filename,
  'help!' => \$help,
  'type=s' => \$type,
  'run!' => \$run,
  'check_all!' => \$check_all,
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





my $fa = $db->get_FileAdaptor;

my $files;
if($filename){
  my $file = $fa->fetch_by_name($filename);
  push(@$files, $file) if($file);
  print STDERR "Failed to find a file using ".$filename."\n" unless($file);
}elsif($type){
  $files = $fa->fetch_by_type($type);
}else{
  $files = $fa->fetch_all;
}



foreach my $file(@$files){
  next if($file->md5 && !$check_all);
  next if($filename && $file->name ne $filename);
  if($run){
    my $md5 = run_md5($file->name);
    next if($file->md5 && ($md5 eq $file->md5));
    $file->md5($md5);
    my $history = ReseqTrack::History->new(
      -other_id => $file->dbID,
      -table_name => 'file',
      -comment => "Updating md5 value",
        );
    $file->history($history);
    $fa->update($file);
  }
}

sub useage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/file/run_md5.pl

=head1 SYNOPSIS

This script will run md5sum to update the md5 value of the given files, by default
it will skip any file which already has an md5

=head1 OPTIONS

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the 
    standard port for mysql-g1kdcc.ebi.ac.uk
-file, the full path to the file you wish to update
-type, To reduce database load the script fetches all the objects then updates the
ones it has md5s for. If the type is specified, then only objects with this type
are fetched so the set to be updated must all be of this typ
-run, Binary flag to indicate sql should be executed
-check_all, this forces md5sum run on everything fetched regardless if it already
has an md5 defined
-help, Binary flag to indicate if the help should be printed out

=head1 Examples

This would update all the files in the md5 file regardless of the type

perl ReseqTrack/scripts/file/update_md5.pl -dbhost myhost -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname my_database -md5_file files.md5

will only update those objects which have type FASTQ

perl ReseqTrack/scripts/file/update_md5.pl -dbhost myhost -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname my_database -md5_file files.md5 -type FASTQ

=cut

