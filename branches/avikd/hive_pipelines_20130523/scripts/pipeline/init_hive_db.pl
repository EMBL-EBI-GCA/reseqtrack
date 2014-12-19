#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use ReseqTrack::Tools::GeneralUtils qw(execute_system_command);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use POSIX qw(strftime);
use Cwd qw(abs_path);

my ($dbhost, $dbuser, $dbpass, $dbport, $dbname);
my ($hive_host, $hive_user, $hive_pass, $hive_port, $hive_name);
my $pipeline_name;
my ($ensembl_cvs_dir, $ensembl_hive_version);

my %options;
&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=i'      => \$dbport,
  'hive_host=s'      => \$hive_host,
  'hive_name=s'      => \$hive_name,
  'hive_user=s'      => \$hive_user,
  'hive_pass=s'      => \$hive_pass,
  'hive_port=i'      => \$hive_port,
  'pipeline_name=s'  => \$pipeline_name,
  'ensembl_cvs_dir=s'  => \$ensembl_cvs_dir,
  'ensembl_hive_version=s'  => \$ensembl_hive_version,
  );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $pipeline = $db->get_PipelineAdaptor->fetch_by_name($pipeline_name);
throw("did not find pipeline with name $pipeline_name") if !$pipeline;

my $init_script = "$ensembl_cvs_dir/ensembl-hive/scripts/init_pipeline.pl";
throw("script doesn't exist $init_script") if ! -f $init_script;

$ensembl_hive_version //= (split(/\/+/, abs_path($ensembl_cvs_dir)))[-1];

$hive_name //= join('_', $dbname, $pipeline_name, strftime("%Y%m%d_%H%M", localtime));

my @cmd_words = ('perl', $init_script, $pipeline->config_module);
if ($pipeline->config_options) {
  push(@cmd_words, $pipeline->config_options);
}
push(@cmd_words, '-host', $hive_host);
push(@cmd_words, '-password', $hive_pass) if $hive_pass;
push(@cmd_words, '-pipeline_db', "-user=$hive_user");
push(@cmd_words, '-pipeline_db', "-port=$hive_port");
push(@cmd_words, '-pipeline_db', "-dbname=$hive_name");

push(@cmd_words, '-reseqtrack_db', "-host=$dbhost");
push(@cmd_words, '-reseqtrack_db', "-user=$dbuser");
push(@cmd_words, '-reseqtrack_db', "-port=$dbport");
push(@cmd_words, '-reseqtrack_db', "-dbname=$dbname");
push(@cmd_words, '-reseqtrack_db', "-pass=$dbpass") if $dbpass;

execute_system_command(join(' ', @cmd_words));

my $url = "mysql://$hive_host:$hive_port/$hive_name";

my $hive_db = ReseqTrack::HiveDB->new(
    -url => $url, -pipeline => $pipeline,
    -hive_version => $ensembl_hive_version,
);

$db->get_HiveDBAdaptor->store($hive_db);

=pod

=head1 NAME

reseqtrack/scripts/pipeline/init_hive_db.pl

=head1 SYNOPSIS

This script initialises a hive database using the configuration defined in the pipeline table
  

=head1 OPTIONS

  database options for the reseqtrack database (must exist and have a pipeline defined):

      -dbhost, the name of the mysql-host
      -dbname, the name of the mysql database
      -dbuser, the name of the mysql user
      -dbpass, the database password if appropriate
      -dbport, the port the mysql instance is running on

  database options for the hive database (to be created):

      -hive_host, the name of the mysql-host
      -hive_user, the name of the mysql user
      -hive_pass, the database password if appropriate
      -hive_port, the port the mysql instance is running on
      (optional) -hive_name, the name of the mysql database, defaults to something sensible and unique

  other options:

  -pipeline_name, refers to a name in the pipeline table
  -ensembl_cvs_dir, path to the ensembl api
  -ensembl_hive_version, (optional) default is to infer it from ensembl_cvs_dir

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"
    $HIVE_DB_OPTS="-hive_host mysql-host -hive_user rw_user -hive_pass **** -hive_port 4175"


  perl reseqtrack/process/run_picard.pl $DB_OPTS $HIVE_DB_OPTS
    -ensembl_cvs_dir /path/to/ensembl-100
    -pipeline_name alignment

=cut
