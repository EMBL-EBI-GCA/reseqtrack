#!/sw/bin/perl -w

use strict;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileSystemUtils qw( get_lines_from_file );
use Getopt::Long;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;

my $file;
my $help;

&GetOptions(
    'dbhost=s' => \$dbhost,
    'dbname=s' => \$dbname,
    'dbuser=s' => \$dbuser,
    'dbpass=s' => \$dbpass,
    'dbport=s' => \$dbport,
    'file=s'   => \$file,
    'help!'    => \$help,
);

if ($help) {
    perldocs();
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);

my $sida = $db->get_StudyIDAdaptor;
my $existing_study_ids = $sida->fetch_all;
my %existing_study_ids;
foreach my $study_id (@$existing_study_ids) {
    $existing_study_ids{$study_id} = 1;
}


my $lines = get_lines_from_file($file);
my %study_ids;
foreach my $line (@$lines) {
    next $line if ($line =~ /^#/);
    my @study_ids = split(/\s+/, $line);

    ID:
    foreach my $study_id (@study_ids) {
        next ID if $existing_study_ids{$study_id};
        $sida->store($study_id);
        $existing_study_ids{$study_id} = 1;
    }
}



sub perldocs {
    exec('perldoc', $0);
    exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/files/study_ids.pl

=head1 SYNOPSIS

    This script takes a config file and parses it to load the
    study_id table

=head2 CONFIG FILE FORMAT

    Lines beginning with '#' are ignored
    Lines may contain multiple study_ids (strings)
    Each study_id is separated by whitespace

    e.g.
        # some G1K study_ids
        SRP004068 SRP000805 SRP004055 SRP001515 SRP004062 SRP000547
        SRP004078 SRP000546 SRP004364 SRP001517 SRP004064 SRP001293

=head2 OPTIONS

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on
        -file, name of file to read
        -help, flag to print this help and exit


=head1 Example:


$DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"

    This will load the study_ids specifed in study_ids.conf:
    perl ReseqTrack/scripts/run_meta_info/load_study_ids.pl  $DB_OPTS -file /path/to/study_ids.conf

=cut

