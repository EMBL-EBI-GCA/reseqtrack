#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use ReseqTrack::Tools::GeneralUtils qw(execute_system_command);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::RunHive;
use Cwd qw(getcwd);
use DateTime::Format::MySQL qw(parse_datetime);

my ($dbhost, $dbuser, $dbpass, $dbport, $dbname);
my ($hive_user, $hive_pass);
my ($pipeline_name, $hive_db_id, $reseed, $loop, $run, $sync, $hive_log_dir);
my $ensembl_cvs_dir;

my %options;
&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=i'      => \$dbport,
  'hive_user=s'      => \$hive_user,
  'hive_pass=s'      => \$hive_pass,
  'pipeline_name=s'  => \$pipeline_name,
  'hive_db_id=s'  => \$hive_db_id,
  'reseed!'  => \$reseed,
  'ensembl_cvs_dir=s'  => \$ensembl_cvs_dir,
  'loop!' => \$loop,
  'run!' => \$run,
  'hive_log_dir=s'  => \$hive_log_dir,
  );
$reseed //= 1;
$loop //= 1;
$run //= 1;
$sync //= 0;

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

throw("must specify a hive_db_id or a pipeline_name")
    if ! defined $hive_db_id && ! defined $pipeline_name;
my $hive_db;
if (defined $hive_db_id) {
  $hive_db = $db->get_HiveDBAdaptor->fetch_by_dbID($hive_db_id);
  throw("did not find a hive_db with id $hive_db_id") if !$hive_db;
  if (defined $pipeline_name) {
    throw("inconsistent pipeline names: $pipeline_name ".$hive_db->pipeline->name)
      if $pipeline_name ne $hive_db->pipeline->name;
  }
}
else {
  my $pipeline = $db->get_PipelineAdaptor->fetch_by_name($pipeline_name);
  throw("did not find pipeline with name $pipeline_name") if !$pipeline;

  my $hive_dbs = $db->get_HiveDBAdaptor->fetch_by_pipeline($pipeline, 1);
  throw("did not find a hive_db for $pipeline_name") if !@$hive_dbs;

  # preferentially take a pipeline that is already seeded:
  my @seeded_hive_dbs = grep {$_->is_seeded} @$hive_dbs;
  if (@seeded_hive_dbs) {
    $hive_dbs = \@seeded_hive_dbs;
  }

  # preferentially take the newest hive_db:
  ($hive_db) = sort {DateTime::Format::MySQL->parse_datetime($b->created)->epoch <=> DateTime::Format::MySQL->parse_datetime($a->created)->epoch} @$hive_dbs;
}

my $hive_dbname = $hive_db->name;
my $hive_port = $hive_db->port;
my $hive_host = $hive_db->host;
#$url =~ s{^(\w*)://(?:[^/\@]*\@)?}{$1://$hive_user:$hive_pass\@};
print "selected hive_db name=$hive_dbname port=$hive_port host=$hive_host\n";

my $run_hive = ReseqTrack::Tools::RunHive->new(
    -hive_dbname => $hive_dbname,
    -hive_scripts_dir => $ensembl_cvs_dir . '/ensembl-hive/scripts',
    -hive_user => $hive_user, -hive_password => $hive_pass,
    -hive_host => $hive_host, -hive_port => $hive_port,
    );

if ($reseed) {
  if ($hive_db->is_seeded) {
    print "no need to reseed hive_db; it is already seeded\n";
  }
  else {
    $run_hive->run('seed', 'get_seeds', '{"seed_time"=>'.time.'}');

    $hive_db->is_seeded(1);
    $db->get_HiveDBAdaptor->update($hive_db);
  }
}

if ($sync) {
  $run_hive->run('sync');
}
if ($run || $loop) {
  $run_hive->options('loop', $loop);
  if ($hive_log_dir) {
    $run_hive->output_dir($hive_log_dir);
    $run_hive->options('use_log_dir', 1);
  }
  $run_hive->run('run');

}


=pod

=head1 NAME

reseqtrack/scripts/pipeline/run_pipeline.pl

=head1 SYNOPSIS

This script is used to do two things:
    1. create a job in a hive database in order to seed the pipeline (make it ready to run)
    2. run the hive pipeline

=head1 OPTIONS

  database options for the reseqtrack database (i.e. NOT the hive database)

      -dbhost, the name of the mysql-host
      -dbname, the name of the mysql database
      -dbuser, the name of the mysql user
      -dbpass, the database password if appropriate
      -dbport, the port the mysql instance is running on

  database options for the hive database:

      -hive_user, the name of the mysql user
      -hive_pass, the database password if appropriate
      (all other hive_db params are got from the reseqtrack hive_db table)

  other options:

  -pipeline_name, refers to a name in the pipeline table. Must specify either this or -hive_db_id
  -hive_db_id, refers to a dbID in the hive_db table. Must specify either this or -pipeline_name
  -ensembl_cvs_dir, path to the ensembl api
  -reseed, boolean flag to add a new seed job to the pipeline. Default is 1. Use -noreseed to disable.
  -run, boolean flag to run the pipeline. Default is 1. Use -norun to disable.
  -sync, boolean flag to sync the pipeline. Default is 0.
  -loop, boolean flag to run the pipeline continuously until nothing let to be done.  Default is 1. Use -noloop to disable.
  -hive_log_dir, optional path to a directory.
        All stdout and stderr for all hive jobs will be written to this directory
        Unnecessary except for debugging because important messages are stored in the hive db.

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/pipeline/run_pipeline.pl $DB_OPTS $HIVE_DB_OPTS
    -pipeline_name alignment
    -ensembl_cvs_dir /path/to/api

=cut
