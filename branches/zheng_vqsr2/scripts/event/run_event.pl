#!/sw/arch/bin/perl -w

use strict;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::PipelineUtils;
use ReseqTrack::Tools::Exception;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $event_name;
my $input_string;
my $verbose;
my $help;

&GetOptions(
    'dbhost=s'       => \$dbhost,
    'dbname=s'       => \$dbname,
    'dbuser=s'       => \$dbuser,
    'dbpass=s'       => \$dbpass,
    'dbport=s'       => \$dbport,
    'event_name=s'   => \$event_name,
    'input_string=s' => \$input_string,
    'help!'          => \$help,
    'verbose!'       => \$verbose,
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

my $ea = $db->get_EventAdaptor;

my $event = $ea->fetch_by_name($event_name);
print "Have event ".$event."\n";
if(!$event){
throw("Failed to fetch an event using ".$event_name);
}
if (!$input_string) {
    $input_string = get_random_input_string($db, $event);
}
if(!$input_string){
  throw("FAiled to get input string for ".$event->name." ".$event->table_name." ".$event->type);
}
my $cmd_line = create_event_commandline($event, $input_string);
print $cmd_line. "\n";
eval {
    my $exit = system($cmd_line);
    throw($cmd_line . " returned a non zero exit code " . $exit)
      unless ($exit == 0);
};

if ($@) {
    throw("Problem running " . $event->name . " with " . $input_string . " $@");
}

sub perldocs {
    exec('perldoc', $0);
    exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/run_event.pl

=head1 SYNOPSIS

This script will run the commandline for one event on either a given, or a random
input 

=head1 OPTIONS

-dbhost, host name for database instance

-dbname, name of database on instance

-dbuser, name of user

-dbpass, password string if required

-dbport, the port number of the database instance

-event_name, the name of the event in the event table

-input_string, the input to be given to the event, if this isn't specified a random string from the database is selected

-help, binary flag to print out the perl docs

-verbose, binary flag to print out additional statements about the process

=head1 Examples

=head2 Other useful scripts

ReseqTrack/event/run_event_pipeline.pl, this will run all the inputs associated with 
the events in your database using LSF on a farm

=cut

