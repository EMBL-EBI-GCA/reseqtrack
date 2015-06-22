#!/usr/bin/env perl

use strict;
use warnings;
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
my $run;
my $help;
my $type;
my $trim;
my $root;

&GetOptions(
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'trim!' => \$trim,
  'root=s' => \$root,
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



my $fa = $db->get_FileAdaptor;

my $files;
if($type){
  $files = $fa->fetch_by_type($type);
}else{
  $files = $fa->fetch_all;
}

foreach my $file(@$files){
  my $name = $file->name;
  if($trim){
    $name =~ s/$root//;
  }
  print $name."\n";
}

sub useage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

    ReseqTrack/scripts/file/list_files.pl

=head1 SYNOPSIS

    This script will provide a list of all entries in the "file" table. Can specify file type

=head1 OPTIONS

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on, this defaults to 4197 the 
      standard port for mysql-g1kdcc.ebi.ac.uk
    -type, Normally all files are fetched. If the type is specified, then only objects 
      with this type dumped
    -help, Binary flag to indicate if the help should be printed out
    -trim, Binary flag to indicate part of the path should be trimmed
    -root, string, pattern to remove

=head1 Examples

    perl ReseqTrack/scripts/file/list_files.pl -dbhost myhost -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname my_database -type FILTERED_FASTQ
    perl ReseqTrack/scripts/file/list_files.pl -dbhost myhost -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname my_database -trim -root /nfs/1000g-archive/vol1/ftp/


=cut

