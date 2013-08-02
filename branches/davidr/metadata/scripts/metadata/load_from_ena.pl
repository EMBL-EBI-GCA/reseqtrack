#! /usr/bin/env perl
use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::ERAUtils qw(get_erapro_conn);

use ReseqTrack::Tools::Metadata::EnaReadInfo;
use ReseqTrack::Tools::GeneralUtils qw(current_date);

use ReseqTrack::Tools::UpdateMetaData;

my %db_params;
my @era_params;
my $load_new;
my $update_existing;
my $help;
my $verbose;
my $force_update;
my $seed_from_old_study_table;
my @manipulator_modules;
my @studies_to_add;
my $clob_read_length = 66000;
my $skip_run_stats;
my @target_studies;
my @target_types;
my $log_dir;

&GetOptions(
  'dbhost=s'                  => \$db_params{-host},
  'dbname=s'                  => \$db_params{-dbname},
  'dbuser=s'                  => \$db_params{-user},
  'dbpass=s'                  => \$db_params{-pass},
  'dbport=s'                  => \$db_params{-port},
  'era_dbuser=s'              => \$era_params[0],
  'era_dbpass=s'              => \$era_params[1],
  'new_study=s'               => \@studies_to_add,
  'load_new!'                 => \$load_new,
  'update_existing!'          => \$update_existing,
  'help!'                     => \$help,
  'verbose!'                  => \$verbose,
  'force_update!'             => \$force_update,
  'manipulator_module=s'      => \@manipulator_modules,
  'clob_read_length=i'        => \$clob_read_length,
  'skip_run_stats!'           => \$skip_run_stats,
  'study=s'                   => \@target_studies,
  'type=s'                    => \@target_types,
  'seed_from_old_study_table' => \$seed_from_old_study_table,
  'log_dir=s'                 => \$log_dir,
);

my $era_db = get_erapro_conn(@era_params);
$era_db->dbc->db_handle->{LongReadLen} = $clob_read_length;

my $reseq_db = ReseqTrack::DBSQL::DBAdaptor->new(%db_params);
$reseq_db->dbc->db_handle->{'AutoCommit'} = 0;

if ($seed_from_old_study_table) {
  my $studyIDAdaptor = $reseq_db->get_StudyIDAdaptor();
  my $study_ids      = $studyIDAdaptor->fetch_all();
  push @studies_to_add, @$study_ids;
}

my @manipulators;
for my $module (@manipulator_modules) {
  my $file = "$module.pm";
  $file =~ s{::}{/}g;
  eval { require "$file"; };
  if ($@) {
    throw("cannot load $file: $@");
  }

  push @manipulators,
    $module->new( -era_db => $era_db, -reseq_db => $reseq_db );
}

my $log_fh;
if ($log_dir) {
  my $ident         = int( rand(10000) );
  my $date          = current_date();
  my $tmp_name      = "update_runmetainfo." . $date . "." . $ident . ".$$.log";
  my $log_file_path = $log_dir . "/" . $tmp_name;
  open( $log_fh, '>', $log_file_path )
    or die("Could not write to $log_file_path: $!");
}

my $updater = ReseqTrack::Tools::UpdateMetaData->new(
  -dcc_db          => $reseq_db,
  -era_db          => $era_db,
  -verbose         => $verbose,
  -log_fh          => $log_fh,
  -manipulators    => \@manipulators,
  -types           => \@target_types,
  -target_studies  => \@target_studies,
  -use_default_rsm => !$skip_run_stats,
);

$updater->load_new_studies(@studies_to_add) if (@studies_to_add);
$updater->update_from_era( $load_new, $update_existing, $force_update );

$reseq_db->dbc->db_handle->commit();

