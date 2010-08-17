#!/usr/local/bin/perl -w

=pod

=head1 NAME

ReseqTrack/scripts/event/runner.pl

=head1 SYNOPSIS

This script is used by the event pipeline system to run the event associated with
a particular job or set of jobs

=head1 OPTIONS

-dbhost, host name for database instance

-dbname, name of database on instance

-dbuser, name of user

-dbpass, password string if required

-dbport, the port number of the database instance

=head1 Examples

perl /homes/laura/code/1000genomes/ReseqTrack/scripts/event/runner.pl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -dbname lec_era_meta_info -name fetch
_archive_fastq 15258

=head2 Other useful scripts

ReseqTrack/event/run_event_pipeline.pl, this is the script for creating job objects to run on the farm and running them.

=cut

use strict;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::EventComplete;
use ReseqTrack::Tools::Exception;
use Getopt::Long;
use Sys::Hostname;
use File::Basename;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 3306;
my $dbname;
my $logic_name;
my @input_strings;
my $array;
print "Commandline runner.pl " . join(" ", @ARGV) . "\n";
&GetOptions(
    'dbhost=s' => \$dbhost,
    'dbname=s' => \$dbname,
    'dbuser=s' => \$dbuser,
    'dbpass=s' => \$dbpass,
    'dbport=s' => \$dbport,
    'name=s'   => \$logic_name,
);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);

throw("runner.pl can't run with a database connection") if (!$db);
throw("runner.pl can't run withtout a logic_name")      if (!$logic_name);
throw("runner.pl needs an input string, this is not a jobarray job")
  if (!$array && (!@ARGV || @ARGV == 0));

my $array_index;
my $job;
my $ja = $db->get_JobAdaptor;
my $aa = $db->get_EventAdaptor;
my @jobs;
print "Have ARGV array @ARGV\n";
foreach my $dbID (@ARGV) {
    $job = $ja->fetch_by_dbID($dbID);
    throw("Failed to fetch job with " . $dbID) if (!$job);
    my $hostname = [split(/\./, hostname())];
    my $host = shift(@$hostname);
    $job->host($host);
    my $id = $ENV{'LSB_JOBID'};
    unless($id){
      job_failed(
            $job,
            "Can't run " 
              . $job . " "
              . $job->dbID . " "
              . $job->input_string
              . " without an submission id",
            $ja
        );
      throw("Failed to get lsf job id");
    }
    $job->submission_id($id);
    unless($job->submission_id){
      job_failed(
            $job,
            "Can't run " 
              . $job . " "
              . $job->dbID . " "
              . $job->input_string
              . " without an submission id",
            $ja
        );
      throw("Don't have a submission id associated with ".$job->dbID);
    }
    $ja->update($job);
    $job->current_status('WAITING');
    $ja->set_status($job);
    $ja->update($job);
    push(@jobs, $job);
}
print "Setting disconnect when inactive to 1\n";
$db->dbc->disconnect_when_inactive(1);
JOB: foreach my $job (@jobs) {
    my $input_string = $job->input_string;
    $job->current_status('RUNNING');
    $ja->set_status($job);
    my $analysis = $aa->fetch_by_name($logic_name);
    if (!$analysis) {
        job_failed(
            $job,
            "Can't run " 
              . $job . " "
              . $job->dbID . " "
              . $job->input_string
              . " without an analysis",
            $ja
        );
    }
    print "Have analysis " . $analysis . "\n";
    my $cmd = $analysis->program . " ";
    $cmd .= $analysis->options . " "    if ($analysis->options);
    $cmd .= $analysis->input_flag . " " if ($analysis->input_flag);
    $cmd .= $input_string;
    print $cmd. "\n";
    my $exit;
    eval {
        print "\n*****" . $cmd . " output******\n";
        $exit = system($cmd);
        if ($exit >= 1) {
            throw($cmd . " returned exit code " . $exit);
        }
        print "**********\n";
    };
    if ($@) {
        my $warning = $cmd . " failed with $@ error";
        job_failed($job, $warning, $ja);
    } elsif ($exit != 0) {
        my $warning = $cmd . " failed with " . $exit . " exit code";
        job_failed($job, $warning, $ja);
        next JOB;
    } else {
        $job->current_status('SUCCESSFUL');
        $ja->set_status($job);
    }
    my $ca               = $db->get_EventCompleteAdaptor;
    my $other_name       = $job->input_string;
    my $completed_string = ReseqTrack::EventComplete->new(
        -event      => $analysis,
        -other_name => $job->input_string,
        -success    => 0,
        -adaptor    => $ca,
    );
    eval {
        $ca->store($completed_string)
          if ($job->current_status eq 'SUCCESSFUL');
    };
    if ($@) {
        my $warning =
          "Failed to store " . $input_string . " in completed_string $@";
        job_failed($job, $warning, $ja, 'FAIL_NO_RETRY');
    }

    if ($job->current_status eq 'SUCCESSFUL') {
        eval { $ja->remove($job); };
        throw("Failed to remove " . $job . " from " . $dbname . " $@")
          if ($@);
    }
}

sub job_failed {
    my ($job, $warning, $ja, $status) = @_;
    $status = 'FAILED' if (!$status);
    print STDERR $warning . "\n";
    $job->current_status($status);
    $ja->set_status($job);
    $ja->unset_submission_id();
}
