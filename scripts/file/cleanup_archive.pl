#!/usr/local/bin/perl -w

use strict;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::ArchiveUtils;
use ReseqTrack::Tools::Exception;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $loops = 1;
my $sleep = 120;
my $verbose = 0;
my $help ;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'loops=s' => \$loops,
  'sleep=s' => \$sleep,
  'verbose!' => \$verbose,
   'help!' => \$help,
    );

if($help){
  useage();
}

#connect to db
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );


my $aa = $db->get_ArchiveAdaptor;



my $count = 0;
while($count <= $loops){
  my $archives = $aa->fetch_all;
  cleanup_archive($archives, $db, $verbose);
  $count++;
  sleep($sleep) if($count <= $loops);
}

$aa->delete_archive_lock;

sub useage{
  exec('perldoc', $0);
  exit(0);
}


=pod

=head1 NAME

  ReseqTrack/scripts/file/cleanup_archive.pl

=head1 SYNOPSIS

    This script silently monitors the progress of file archiving. If a file has 
    been successfully moved the entry in the "archive" table relating to that 
    file is then removed. 

=head1 OPTIONS

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on, this defaults to 4197 the 
     standard port for mysql-g1kdcc.ebi.ac.uk
    -sleep, the amount of time to wait  between cleaning loops ( default 120 seconds)
    -help, Binary flag to indicate if the help should be printed out



=head1 Examples

    perl ReseqTrack/scripts/file/cleanup_archive.pl -dbhost myhost -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname my_database 

=cut
