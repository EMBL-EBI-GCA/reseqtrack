#!/usr/local/bin/perl -w

use strict;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::PipelineUtils;
use Getopt::Long;
use File::Copy;
use File::Path;
use Benchmark;
$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;
my $current;
my $finished;
my $current_summary;
my $input;
my $finished_percent;
my $help;
my $event_summary;

&GetOptions(
    'dbhost=s'          => \$dbhost,
    'dbname=s'          => \$dbname,
    'dbuser=s'          => \$dbuser,
    'dbpass=s'          => \$dbpass,
    'dbport=s'          => \$dbport,
    'current!'          => \$current,
    'finished!'         => \$finished,
    'current_summary!'  => \$current_summary,
    'input!'            => \$input,
    'finished_percent!' => \$finished_percent,
    'help!'             => \$help,
    'event_summary!'    => \$event_summary,
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

my $ja = $db->get_JobAdaptor;
my $ca = $db->get_EventCompleteAdaptor;
my $aa = $db->get_EventAdaptor;

my ($jobs, $analyses, $input_hash);

$jobs       = $ja->fetch_all if ($current || $current_summary);
$analyses   = $aa->fetch_all;
$input_hash = get_input_counts($analyses, $db)
  if ($input || $finished_percent);

my %status_count;
my %logic_status_count;
my $completed_event = get_completed_event_hash($db);
my %input_string;
my %analyses;
my %event_hash;
foreach my $analysis (@$analyses) {
    $event_hash{$analysis->name} = $analysis;
}

foreach my $job (@$jobs) {
    if (!$status_count{$job->current_status}) {
        $status_count{$job->current_status} = 1;
    } else {
        $status_count{$job->current_status}++;
    }
    if (!$logic_status_count{$job->event->name}) {
        $logic_status_count{$job->event->name} = {};
    }
    if (!$logic_status_count{$job->event->name}{$job->current_status}) {
        $logic_status_count{$job->event->name}{$job->current_status} = 1;
    } else {
        $logic_status_count{$job->event->name}{$job->current_status}++;
    }
}

if ($input) {
    print "\n#INPUT\n";
    my %printed;
    foreach my $event_name (keys(%$input_hash)) {
        my $event = $event_hash{$event_name};
        next if ($printed{$event_name});
        print $event->type . " "
          . $event->table_name . " "
          . $input_hash->{$event_name} . "\n";
        $printed{$event_name} = 1;
    }
    print "\n";
}

if ($current) {
    print "\n#CURRENT\n";
    foreach my $name (keys(%logic_status_count)) {
        foreach my $status (keys(%{$logic_status_count{$name}})) {
            print $name. " " 
              . $status . " "
              . $logic_status_count{$name}->{$status} . "\n";
        }
    }
    print "\n";
}
if ($current_summary) {
    print "\n#CURRENT SUMMARY\n";
    foreach my $status (keys(%status_count)) {
        print $status. " " . $status_count{$status} . "\n";
    }
    print "\n";
}
if ($finished) {
    print "\n#COMPLETED\n";
    foreach my $logic (keys(%$completed_event)) {
        my $event = $event_hash{$logic};
        my $count = $completed_event->{$logic}->{$event->type};
        print $logic. "\t" . $event->type . "\t" . $count . "\n";
    }
    print "\n";
}
if ($finished_percent) {
    print "\n#PERCENT COMPLETE\n";
    foreach my $logic (keys(%$completed_event)) {
        my $event = $event_hash{$logic};
        throw("Don't have an event for " . $logic) unless ($event);
        my $num_complete    = $completed_event->{$logic}->{$event->type};
        my $number_input    = $input_hash->{$event->name};
        my $percent         = ($num_complete / $number_input) * 100 
	  if($num_complete && $number_input);
	$percent = 0 unless($percent);
        my $rounded_percent = sprintf("%.2f", $percent);
        print $logic. " " . $rounded_percent . "%\n";
    }
    print "\n";
}
if ($event_summary) {
    print "\nEvent Summary\n";
    foreach my $event (@$analyses) {
        print $event->name . " "
          . $event->farm_options . " "
          . $event->table_name . " "
          . $event->type . "\n";
    }
}

sub get_input_counts {
    my ($events, $db) = @_;
    my %event_count_hash;
    foreach my $event (@$events) {
        if ($event->table_name eq 'file') {
            my $sql = "select count(*) as count, type from file group by type";
            my $sth = $db->dbc->prepare($sql);
            $sth->execute;
            while (my ($count, $type) = $sth->fetchrow) {
                $event_count_hash{$event->name} = $count
                  if ($event->type eq $type);
            }
        } elsif ($event->table_name eq 'collection') {
            my $sql = "select count(*), type from collection group by type";
            my $sth = $db->dbc->prepare($sql);
            $sth->execute;
            while (my ($count, $type) = $sth->fetchrow) {
                $event_count_hash{$event->name} = $count
                  if ($event->type eq $type);
            }
        } elsif ($event->table_name eq 'run_meta_info') {
            my $sql = "select count(*) from run_meta_info";
            my $sth = $db->dbc->prepare($sql);
            $sth->execute;
            while (my ($count) = $sth->fetchrow) {
                $event_count_hash{$event->name} = $count;
            }
        } elsif ($event->table_name eq 'input_string') {
            my $sql = "select type, count(*) from input_string group by type";
            my $sth = $db->dbc->prepare($sql);
            $sth->execute;
            while (my ($type, $count) = $sth->fetchrow) {
                $event_count_hash{$event->name} = $count
                  if ($event->type eq $type);
            }
        }
    }
    return \%event_count_hash;
}

sub get_completed_event_hash {
    my ($db) = @_;
    my $sql =
"select event.name, event.type, count(distinct(other_id)) from event, event_complete where event.event_id = event_complete.event_id group by event.name , event.type";
    my $sth = $db->dbc->prepare($sql);
    $sth->execute;
    my %hash;
    while (my ($name, $type, $count) = $sth->fetchrow) {
        $hash{$name}->{$type} = $count;
    }
    return \%hash;
}

sub perldocs {
    exec('perldoc', $0);
    exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/event/monitor_event_pipeline.pl

=head1 SYNOPSIS

This script provides a summary of the input ids, existing jobs and sucessfully 
completed jobs in the system

=head1 OPTIONS

-dbhost, host name for database instance

-dbname, name of database on instance

-dbuser, name of user

-dbpass, password string if required

-dbport, the port number of the database instance

-current, a summary of the jobs currently in the system by event name and status

-finished, a summary of completed jobs from the event complete table, grouped by event_name

-current_summary, a summary of jobs currently in the system grouped only by status

-finished_percent, a summary of how many jobs have finish giving percentages of input rather than raw counts

-input, a summary of the number of inputs for each table_name and type

-event_summary, a summary of the events in the system

-help, a binary flag to print out the perl docs

=head1 Examples

 perl ReseqTrack/scripts/event/monitor_event_pipeline.pl  -dbhost mysql-host -dbuser ro_user -dbport 4197 -dbname my_database  -current -finished -input -finished_p

This will produce output like

INPUT
run_meta_info RUNMETAINFO 15461
collection ARCHIVE_FASTQ 13707


CURRENT
archive_fastq_sanity FAILED 99
fetch_archive_fastq FAILED 1567
fetch_archive_fastq SUBMITTED 1


COMPLETED
archive_fastq_sanity 13608
fetch_archive_fastq 13893


PERCENT COMPLETE
archive_fastq_sanity 99%
fetch_archive_fastq 90%

=head2 Other useful scripts

ReseqTrack/event/run_event_pipeline.pl is a script which can submit jobs to the farm and work out what event/input combinations can be run

=cut

