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

&GetOptions(
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'md5_file=s' => \$md5_file,
  'help!' => \$help,
  'type=s' => \$type,
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

my $fh;
if($md5_file){
  open(FH, ">".$md5_file) or throw("FAiled to open ".$md5_file);
  $fh = \*FH;
}else{
  $fh = \*STDOUT;
}

my $fa = $db->get_FileAdaptor;

my $files;
if($type){
  $files = $fa->fetch_by_type($type);
}else{
  $files = $fa->fetch_all;
}

foreach my $file(@$files){
  print $fh $file->md5."\t".$file->name."\n";
}

sub useage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/file/dump_md5.pl

=head1 SYNOPSIS

This script will dump an md5 file

=head1 OPTIONS

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the 
    standard port for mysql-g1kdcc.ebi.ac.uk
-md5_file, the file to write the md5s to
-type, To reduce database load the script fetches all the objects then updates the
ones it has md5s for. If the type is specified, then only objects with this type
are fetched so the set to be updated must all be of this type
-help, Binary flag to indicate if the help should be printed out

=head1 Examples

perl ReseqTrack/scripts/file/dump_md5.pl -dbhost myhost -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname my_database -md5_file files.md5 -type FASTQ

dump the md5s of all files with type FASTQ into afile named files.md5

=cut

