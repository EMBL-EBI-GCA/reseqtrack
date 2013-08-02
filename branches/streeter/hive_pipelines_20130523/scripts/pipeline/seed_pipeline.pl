#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use ReseqTrack::Tools::GeneralUtils qw(execute_system_command);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;

my ($dbhost, $dbuser, $dbpass, $dbport, $dbname);
my ($hive_user, $hive_pass);
my $pipeline_name;
my $ensembl_cvs_dir;
my $run_pipeline;

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
  'ensembl_cvs_dir=s'  => \$ensembl_cvs_dir,
  'run_pipeline!' => \$run_pipeline,
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

my $seed_script = "$ensembl_cvs_dir/ensembl-hive/scripts/seed_pipeline.pl";
throw("script doesn't exist $seed_script") if ! -f $seed_script;

my $hive_dbs = $db->get_HiveDBAdaptor->fetch_by_pipeline($pipeline, 1);
throw("did not find a hive_db for $pipeline_name") if !@$hive_dbs;
throw("a pipeline is already seeded for $pipeline_name")
  if scalar grep {$_->is_seeded} @$hive_dbs;
my $hive_db = $hive_dbs->[0];
my $url = $hive_db->url;
$url =~ s{^(\w*)://(?:[^/\@]*\@)?}{$1://$hive_user:$hive_pass\@};


my @cmd_words = ('perl', $seed_script);
push(@cmd_words, '-url', q(').$url.q('));
push(@cmd_words, '-logic_name', 'get_seeds');
push(@cmd_words, '-input_id', q('{"seed_time"=>).time.q(}'));

execute_system_command(join(' ', @cmd_words));

$hive_db->is_seeded(1);
$db->get_HiveDBAdaptor->update($hive_db);
