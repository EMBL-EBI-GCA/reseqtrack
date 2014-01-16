#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileUtils qw( assign_type calculate_comment);
use ReseqTrack::Tools::FileSystemUtils qw( get_lines_from_file);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::File;
use ReseqTrack::History;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $list_file;
my $default_host;
my $update_existing = 0;
my $help;
&GetOptions(
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'list_file=s' => \$list_file,
  'host=s' => \$default_host,
  'update_existing!' => \$update_existing,
  'help!' => \$help,
   );

my $type_placeholder = 'unassigned';

if($help){
  useage();
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);

my %hosts;
foreach my $host (@{$db->get_HostAdaptor->fetch_all}) {
  $hosts{$host->name} = $host;
}

my @files;
my $fa = $db->get_FileAdaptor;

LINE:
foreach my $line (@{get_lines_from_file($list_file)}) {
  next LINE if $line =~ /^#/;
  next LINE if $line !~ /\S/;
  chomp $line;
  my ($path, $size, $md5, $host) = split("\t", $line);
  throw("invalid line: $line") if !$path || ! defined $size || !$md5;
  $host //= $default_host;
  throw("no host for $path") if !$host;
  throw("did not recognise host: $host") if ! defined $hosts{$host};

  my $basename = fileparse($path);
  my $exists = $fa->fetch_by_filename($basename);

  throw("file already exists and not set to update: $path") if @$exists && !$update_existing;
  throw("more than one file already exists: $path") if @$exists > 1;

  $size =~ s/[\s,"]//g;
  throw("did not recognise size of $path") if $size !~ /^\d+$/;

  my $file = ReseqTrack::File->new(
        -name => $path, -size => $size,
        -md5 => $md5, -host => $hosts{$host},
        -type => $type_placeholder,
        );

  if (@$exists) {
    my $old_file = $exists->[0];

    # check if a pipeline has already run on that file:
    my $existing_pipeline_seeds = $db->get_PipelineSeedAdaptor->fetch_by_seed($old_file);
    if (@$existing_pipeline_seeds) {
      throw("cannot update $path because pipelines have already run on this file: "
          . join(' ', map {$_->pipeline->name} @$existing_pipeline_seeds));
    }
    foreach my $existing_collection (@{$db->get_CollectionAdaptor->fetch_by_other_id_and_table_name($old_file->dbID, 'file')}) {
      $existing_pipeline_seeds = $db->get_PipelineSeedAdaptor->fetch_by_seed($existing_collection);
      if (@$existing_pipeline_seeds) {
        throw("cannot update $path because pipelines have already run on a collection containing this file: "
            . join(' ', map {$_->pipeline->name} @$existing_pipeline_seeds));
      }
    }

    $file->type($old_file->type);
    my $comment = calculate_comment($file, $old_file);
    if (!$comment) {
      print "skipping $path because no difference between new and old";
      next LINE;
    }
    my $history = ReseqTrack::History->new(
      -other_id => $old_file->dbID,
      -table_name => 'file',
      -comment => $comment,
        );
    $file->dbID($old_file->dbID);
    $file->history($history);
  }

  push (@files, $file);
}

assign_type(\@files, $db);
foreach my $file (@files) {
  throw("no file type assigned for ".$file->name) if $file->type eq $type_placeholder;
}

foreach my $file (@files) {
  $fa->store($file, $update_existing);
}


sub useage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/files/load_file_table.pl

=head1 SYNOPSIS

This script will take a list of file information and load it into
the file table of a ReseqTrack database

This script was designed for loading foreign files as a precursor to the file release pipeline

=head1 OPTIONS

Database options

This is the database the objects will be written to

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the standard
port for mysql-g1kdcc.ebi.ac.uk

File list option

-list_file, a file containing a list of files path details
        Tab delimited, ignores lines starting with '#'
        column 1 = path
        column 2 = size
        column 3 = md5
        column 4 = host name (optional, overrides -host option if it is set)

Other options

-host, This is the host name used if host is not present in the list_file

-update_existing, This means if the store method finds a file with the same name it
will use the update method rather than the store method

=head1 Examples

perl loat_file_table.pl $DB_OPTS -list_file my_list_file

where my_list_file looks like this:

#name	size	md5	host
/path/to/file.ext	12345	e0868bea88651a07a5a9e77e7448994e	partner_1

=cut

