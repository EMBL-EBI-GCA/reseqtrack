#! /usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::ERAUtils qw(get_erapro_conn);

use ReseqTrack::Tools::Metadata::EnaReadInfo;
use ReseqTrack::Tools::GeneralUtils qw(current_date);

use ReseqTrack::Tools::UpdateMetaData;
use ReseqTrack::Tools::Metadata::BaseMetaDataClashCheck;

my %db_params;
my @era_params;
my $load_new;
my $update_existing;
my $help;
my $verbose;
my $force_update;
my $seed_from_old_study_table;
my @add_in_modules;
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
  'add_in=s'                  => \@add_in_modules,
  'clob_read_length=i'        => \$clob_read_length,
  'skip_run_stats!'           => \$skip_run_stats,
  'study=s'                   => \@target_studies,
  'type=s'                    => \@target_types,
  'seed_from_old_study_table' => \$seed_from_old_study_table,
  'log_dir=s'                 => \$log_dir,
);

if ($help) {
  usage();
}

my $era_db = get_erapro_conn(@era_params);
$era_db->dbc->db_handle->{LongReadLen} = $clob_read_length;

my $reseq_db = ReseqTrack::DBSQL::DBAdaptor->new(%db_params);
$reseq_db->dbc->db_handle->{'AutoCommit'} = 0;

if ($seed_from_old_study_table) {
  my $studyIDAdaptor = $reseq_db->get_StudyIDAdaptor();
  my $study_ids      = $studyIDAdaptor->fetch_all();
  push @studies_to_add, @$study_ids;
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
else {
  $log_fh = *STDOUT;
}

my @add_ins;
my $have_clash_check = undef;
for my $module (@add_in_modules) {
  my $file = "$module.pm";
  $file =~ s{::}{/}g;
  eval { require "$file"; };
  if ($@) {
    throw("cannot load $file: $@");
  }

  my $add_in = $module->new(
    -era_db   => $era_db,
    -reseq_db => $reseq_db,
    -log_fh   => $log_fh
  );
  push @add_ins, $add_in;

  if ( $add_in->isa("ReseqTrack::Tools::Metadata::BaseMetaDataClashCheck") ) {
    $have_clash_check = 1;
  }
}
if ( !$have_clash_check ) {
  push @add_ins,
    ReseqTrack::Tools::Metadata::BaseMetaDataClashCheck->new(
    -era_db   => $era_db,
    -reseq_db => $reseq_db,
    -log_fh   => $log_fh
    );
}

my $updater = ReseqTrack::Tools::UpdateMetaData->new(
  -dcc_db          => $reseq_db,
  -era_db          => $era_db,
  -verbose         => $verbose,
  -log_fh          => $log_fh,
  -add_ins         => \@add_ins,
  -target_types    => \@target_types,
  -target_studies  => \@target_studies,
  -use_default_rsm => !$skip_run_stats,
);

$updater->load_new_studies(@studies_to_add) if (@studies_to_add);
$updater->update_from_era( $load_new, $update_existing, $force_update );

$reseq_db->dbc->db_handle->commit();

$updater->report();

sub usage {
  exec( 'perldoc', $0 );
  exit(0);
}

=pod

=head1 NAME

reseqtrack/scripts/metadata/load_from_ena.pl

=head1 SYNOPSIS

This script loads metadata from the ENA database (ERAPRO) into a ReseqTrack database.

=head1 OPTIONS

  database options:

    -dbhost, the name of the mysql-host
    -dbname, the name of the mysql database
    -dbuser, the name of the mysql user
    -dbpass, the database password if appropriate
    -dbport, the port the mysql instance is running on
    -era_dbuser, the name of the user on ERAPRO
    -era_dbpass, the password for ERAPRO

  initial set up options:

	-new_study <study_id>, studies to be added to the ReseqTrack database (can be specified many times)
	-seed_from_old_study_table, flad to populate new_study from the study_id table    

	regular options:
	
	-load_new, load new entries from ERA pro into your ReseqTrack database
	-update_existing, update existing entries in your ReseqTrack database from ERAPRO
	-add_in, package name for a module that will transfrom the information from ERAPRO before recording it in the ReseqTrack database. See ReseqTrack::Tools::Metadata::BaseMetadataAddIn. Can be specified many times.
	-skip_run_stats, data from the ENA has Run statistics. data from the EGA does not, so you set this flag and we won't try to load that information
	-log_dir, a new log file will be created in the specified directory. 
		
	advanced options:
	
	-force_update, update existing entries even when no change is detected (useful if a bug has been fixed, or a add_in module has been changed)
	-clob_read_length <integer>, setting passed to the Oracle driver. You should not need to set this.
	-study <erapro_study_id>, only update given study, rather than the default of all known studies
	-type <type>, only update the given data types, rather than all. valid values are: study sample experiment run. It is possible to cause errors with this if you try to load data without a referenced entity being listed (i.e. runs need experiments and samples, experiments need studies)

	misc options:
	
	-help, display this documentation
	-verbose, print logging information to STDOUT
	-suppress_clash_check, will skip the clash check report, unless you have specifically named an add-in that does this

=head1 Examples

    $DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database -era_user an_era_user -era_pass era_password"

	Transfer studies from an old style ReseqTrack DB to use the new metadata schema (SQL schema must have already been updated)
	 
  perl reseqtrack/metadata/load_from_ena.pl $DB_OPTS -seed_from_old_study_table
	
	Add studies:
	perl reseqtrack/metadata/load_from_ena.pl $DB_OPTS -new_study ERP001466 -new_study ERP001664
	
	Normal invocation for a 1000genomes DB (ENA):
	perl reseqtrack/metadata/load_from_ena.pl $DB_OPTS -load_new -update_existing
		-add_in ReseqTrack::Tools::Metadata::G1KManipulator
		-add_in ReseqTrack::Tools::Metadata::PopulationRulesManipulator

	Normal invocation for a Blueprint DB (EGA):
	perl reseqtrack/metadata/load_from_ena.pl $DB_OPTS -load_new -update_existing
		-skip_run_stats


=cut
