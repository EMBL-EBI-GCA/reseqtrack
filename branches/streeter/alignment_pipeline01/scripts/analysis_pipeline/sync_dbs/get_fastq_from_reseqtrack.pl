#!/sw/arch/bin/perl 

use strict;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use Getopt::Long;
use Data::Dumper;

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
my $type;

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
	    'type=s' => \$type,
	    'help!' => \$help,
	   );


if($help){
    exec('perldoc', $0);
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


my $ha_w = $db_w->get_HostAdaptor;
my %host_ids;
foreach my $host_r (@{$db_r->get_HostAdaptor->fetch_all}) {
  my $host_w = $ha_w->fetch_by_name($host_r->name);
  if ($host_w) {
    $host_ids{$host_r->dbID} = $host_w->dbID;
  }
}

my $rmia_w = $db_w->get_RunMetaInfoAdaptor;
my $ca_r = $db_r->get_CollectionAdaptor;
my $ca_w = $db_w->get_CollectionAdaptor;


RMI:
foreach my $run_meta_info (@{$rmia_w->fetch_all}) {
  my $collection = $ca_r->fetch_by_name_and_type($run_meta_info->run_id, $type);
  next RMI if (!$collection);
  foreach my $file (@{$collection->others}) {
    foreach my $statistics (@{$file->statistics}){
      $statistics->dbID(undef);
    }
    $file->host->dbID($host_ids{$file->host_id});
    $file->dbID(undef);
    $file->empty_history;
  }
  foreach my $statistics(@{$collection->statistics}){
    $statistics->dbID(undef);
  }
  foreach my $history(@{$collection->history}){
    $history->dbID(undef);
  }
  $collection->empty_history;
  $ca_w->store($collection);
}


=pod

=head1 NAME

reseqtrack/scripts/sync_dbs/get_fastq_from_reseqtrack.pl

=head1 SYNOPSIS

This script reads fastq entries from a 'read-only' database and loads them into 'write' database.
Gets fastqs for all entries in the run_meta_info table of the write database.

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

  -type, collection type for the fastq files

  -help, flag to print this help and exit

=head1 Examples

    $DB_R_OPTS="-dbhost_r mysql-host -dbuser_r ro_user -dbpass_r **** -dbport_r 4197 -dbname_r my_database_1"
    $DB_W_OPTS="-dbhost_w mysql-host -dbuser_w rw_user -dbpass_w **** -dbport_w 4197 -dbname_w my_database_2"

  perl reseqtrack/scripts/sync_dbs/get_fastq_from_reseqtrack.pl $DB_R_OPTS $DB_W_OPTS -type=FILTERED_FASTQ

=cut

