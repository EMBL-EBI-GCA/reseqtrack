#!/sw/arch/bin/perl 

use strict;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::RunMetaInfoUtils qw(are_run_meta_infos_identical create_history_for_run_meta_info);
use Getopt::Long;

$| = 1;

my $dbhost_w;
my $dbuser_w;
my $dbpass_w;
my $dbport_w = 4175;
my $dbname_w;
my $dbhost_r;
my $dbuser_r;
my $dbpass_r;
my $dbport_r = 4175;
my $dbname_r;
my $help;
my $store_new;
my $update_existing;
my $log_missing;
my $update_study_ids;
my $all_checks;
my $collection_type;
my $collection_name;
my $verbose;

&GetOptions( 
	    'dbhost_r=s'      => \$dbhost_r,
	    'dbname_r=s'      => \$dbname_r,
	    'dbuser_r=s'      => \$dbuser_r,
	    'dbpass_r=s'      => \$dbpass_r,
	    'dbport_r=s'      => \$dbport_r,
	    'dbhost_w=s'      => \$dbhost_w,
	    'dbname_w=s'      => \$dbname_w,
	    'dbuser_w=s'      => \$dbuser_w,
	    'dbpass_w=s'      => \$dbpass_w,
	    'dbport_w=s'      => \$dbport_w,
	    'collection_type=s' => \$collection_type,
	    'collection_name=s' => \$collection_name,
	    'all_checks!' => \$all_checks,
	    'store_new!' => \$store_new,
	    'update_existing!' => \$update_existing,
	    'log_missing!' => \$log_missing,
	    'update_study_ids!' => \$update_study_ids,
	    'verbose!' => \$verbose,
	    'help!' => \$help,
	   );


if($help){
    exec('perldoc', $0);
    exit(0);
}

if($all_checks){
  $store_new = 1;
  $update_existing = 1;
  $log_missing = 1;
  $update_study_ids = 1;
}

if(!$store_new && !$update_existing && !$log_missing && !$update_study_ids){
  print STDERR "There script has nothing to do; you need to specify at least ".
    "one of -store_new, -update_existing, -log_missing, -update_study_ids ".
      "or -all_checks\n";
  exit(0);
}


my $db_r = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost_r,
  -user   => $dbuser_r,
  -port   => $dbport_r,
  -dbname => $dbname_r,
  -pass   => $dbpass_r,
    );

my $db_w = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost_w,
  -user   => $dbuser_w,
  -port   => $dbport_w,
  -dbname => $dbname_w,
  -pass   => $dbpass_w,
    );

my $sia = $db_w->get_StudyIDAdaptor;

my %study_ids;
if ($update_study_ids) {
  my %study_id_pairs;
  foreach my $existing (@{$sia->fetch_all}) {
    $study_id_pairs{$existing}->{'existing'} = 1;
  }
  foreach my $reference (@{$db_r->get_StudyIDAdaptor->fetch_all}) {
    $study_id_pairs{$reference}->{'reference'} = 1;
  }
  while (my ($study_id, $si_pair) = each %study_id_pairs) {
    my $existing = $si_pair->{'existing'};
    my $reference = $si_pair->{'reference'};
    if ($existing && $reference) {
      $study_ids{$study_id} = 1;
    }
    elsif ($existing && !$reference) {
      print "MISSING study_id $study_id in $dbname_r\n" if ($verbose);
    }
    elsif ($reference && !$existing) {
      print "STORING study_id $study_id in $dbname_w\n" if($verbose);
      $sia->store($study_id);
      $study_ids{$study_id} = 1;
    }
  }
}
else {
  foreach my $study_id (@{$sia->fetch_all}) {
    $study_ids{$study_id} = 1;
  }
}

exit(0) if (!$store_new && !$update_existing && !$log_missing);

my $rmia_w = $db_w->get_RunMetaInfoAdaptor;

my %rmis;
foreach my $study_id (keys %study_ids) {
  foreach my $rmi (@{$rmia_w->fetch_by_study_id($study_id)}) {
    $rmis{$rmi->run_id}->{'existing'} = $rmi;
  }
}

if ($collection_type && $collection_name) {
  my $ca = $db_r->get_CollectionAdaptor;
  my $collection = $ca->fetch_by_name_and_type($collection_name, $collection_type);
  throw("did not find collection with name $collection_name and type $collection_type")
    if (! $collection);
  foreach my $rmi (@{$collection->others}) {
    if ($study_ids{$rmi->study_id}) {
      $rmis{$rmi->run_id}->{'reference'} = $rmi;
    }
  }
}
else {
  my $rmia_r = $db_r->get_RunMetaInfoAdaptor;
  foreach my $study_id (keys %study_ids) {
    foreach my $rmi (@{$rmia_r->fetch_by_study_id($study_id)}) {
      $rmis{$rmi->run_id}->{'reference'} = $rmi;
    }
  }
}

RMI:
while (my ($run_id, $rmi_pair) = each %rmis) {
  my $existing = $rmi_pair->{'existing'};
  my $reference = $rmi_pair->{'reference'};
  if ($update_existing && $existing && $reference) {
    if (! are_run_meta_infos_identical($reference, $existing)) {
      $reference->dbID($existing->dbID);
      my $history = create_history_for_run_meta_info($reference, $existing);
      next RMI if (! $history);
      $reference->history($history);
      print "UPDATING run_id " . $reference->run_id." in $dbname_w: ".$history->comment."\n" if($verbose);
      $rmia_w->update($reference);
    }
  }
  if ($store_new && $reference && !$existing) {
    print "STORING run_id " . $reference->run_id . " in $dbname_w\n" if($verbose);
    $rmia_w->store($reference);
  }
  if (($log_missing || $verbose)  && $existing && !$reference) {
    print "MISSING run_id $run_id in $dbname_r\n";
  }
}

=pod

=head1 NAME

ReseqTrack/scripts/run_meta_info/update_run_meta_info_from_reseqtrack.pl

=head1 SYNOPSIS

This script compares the contents of the run_meta_info table in two different reseqtrack databases

It will summarise any new entries and differences and store them if you
let it.

=head1 OPTIONS

  options for the 'read-only' database:

    -dbhost_r, the name of the mysql-host
    -dbname_r, the name of the mysql database
    -dbuser_r, the name of the mysql user
    -dbpass_r, the database password if appropriate
    -dbport_r, the port the mysql instance is running on

  options for the 'write' database to be updated:

    -dbhost_w, the name of the mysql-host
    -dbname_w, the name of the mysql database
    -dbuser_w, the name of the mysql user
    -dbpass_w, the database password if appropriate
    -dbport_w, the port the mysql instance is running on

  other options:

  -store_new, flag to store new run_meta_info from the read-only database

  -update_existing, flag to update run_meta_info when there is a difference between the
  read-only and write databases

  -log_missing, flag to print a message when run_meta_info from the write database is missing
  from the read-only database

  -update_study_ids, flag to store new study_ids from the read-only database. Study_ids in the
  write database will be ignored if they are missing from the read-only database.

  -all_checks, this switches on all the different checks

  -collection_name, if set then only gets run_meta_info from this collection
  -collection_type, if set then only gets run_meta_info from this collection

  -verbose, flag to print details of the differences between the databases

  -help, flag to print this help and exit

=head1 Examples

    $DB_R_OPTS="-dbhost_r mysql-host -dbuser_r ro_user -dbpass_r **** -dbport_r 4197 -dbname_r my_database_1"
    $DB_W_OPTS="-dbhost_w mysql-host -dbuser_w rw_user -dbpass_w **** -dbport_w 4197 -dbname_w my_database_2"

  perl reseqtrack/scripts/sync_dbs/update_run_meta_info_from_reseqtrack.pl $DB_R_OPTS $DB_W_OPTS -all_checks -collection_name=exome -collection_type=STUDY_TYPE

=cut

