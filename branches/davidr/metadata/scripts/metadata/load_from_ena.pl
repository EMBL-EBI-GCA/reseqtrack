#! /usr/bin/env perl
use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::ERAUtils;

my %db_params;
my @era_params;
my $load_new;
my $update_existing;
my $help;
my $verbose;

&GetOptions(
	'dbhost=s'     => \$db_params{-host},
	'dbname=s'     => \$db_params{-dbname},
	'dbuser=s'     => \$db_params{-user},
	'dbpass=s'     => \$db_params{-pass},
	'dbport=s'     => \$db_params{-port},
	'era_dbuser=s' => \$era_params[0],
	'era_dbpass=s' => \$era_params[1],

	'load_new!'        => \$load_new,
	'update_existing!' => \$update_existing,
	'help!'            => \$help,
	'verbose!'         => \$verbose,
);

my $era_db = get_erapro_conn(@era_params);
my $reseq_db = ReseqTrack::DBSQL::DBAdaptor->new(%db_params);

print Dumper($era_db->get_StudyAdaptor()->fetch_by_source_id('ERP001663'));


#my $run_id = 'ERR248910'
#'ERX223448'
#'ERS208290'