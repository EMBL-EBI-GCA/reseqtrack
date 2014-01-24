#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use ReseqTrack::Tools::GeneralUtils qw(delete_lock_string is_locked create_lock_string);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use DateTime::Format::MySQL;

my ($dbhost, $dbuser, $dbpass, $dbport, $dbname, $force, $futile);
my $pipeline_name;

my %options;
&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=i'      => \$dbport,
  'pipeline_name=s'  => \$pipeline_name,
  'force!' => \$force,
  'futile!' => \$futile,
  );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

throw("the futile flag must be defined if using -force") if $force && ! defined $futile;

my $meta_adaptor = $db->get_MetaAdaptor;
my $lock_string = 'pipeline.lock';
eval{ is_locked($lock_string, $meta_adaptor);};
if ($@) {
  throw("ReseqTrack database is locked with $lock_string in meta table");
}
create_lock_string($lock_string, $meta_adaptor);

my $pipeline = $db->get_PipelineAdaptor->fetch_by_name($pipeline_name);
throw("did not find pipeline with name $pipeline_name") if !$pipeline;

my $hive_dbs = $db->get_HiveDBAdaptor->fetch_by_pipeline($pipeline);
my $ps_a = $db->get_PipelineSeedAdaptor;
HIVEDB:
foreach my $hive_db (grep {!$_->retired} @$hive_dbs) {
  my $pipeline_seeds = $db->get_PipelineSeedAdaptor->fetch_running_by_hive_db($hive_db);
  if (@$pipeline_seeds) {
    if (!$force) {
      warn("will not retire ".$hive_db->name." because pipeline seeds are still marked as running and -force is not set");
      next HIVEDB;
    }
    foreach my $ps (@$pipeline_seeds) {
      $ps_a->update_failed($ps, $futile);
    }
  }

  $db->get_HiveDBAdaptor->retire($hive_db);
}

delete_lock_string($lock_string, $meta_adaptor);


=pod

=head1 NAME

reseqtrack/scripts/pipeline/retire_hive_db.pl

=head1 SYNOPSIS

This script marks an old hive database as 'retired' in the hive_db table of a reseqtrack database.
This is so the hive database will not be used any more for running pipelines
  

=head1 OPTIONS

  database options for the reseqtrack database (i.e. NOT the hive database)

      -dbhost, the name of the mysql-host
      -dbname, the name of the mysql database
      -dbuser, the name of the mysql user
      -dbpass, the database password if appropriate
      -dbport, the port the mysql instance is running on

  other options:

  -pipeline_name, refers to a name in the pipeline table
  -force, binary flag. This is needed if the hive_db is still marked as running some pipeline seeds.  These pipeline seeds will be marked as failed.
  -futile, binary flag.  For use with the -force flag.  Failed pipeline seeds can also be marked futile for that pipeline so they will not be retried by another hive_db.

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

  perl reseqtrack/pipeline/retire_hive_db.pl $DB_OPTS
    -pipeline_name alignment
    -force -nofutile

=cut

