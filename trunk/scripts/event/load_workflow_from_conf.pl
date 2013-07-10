#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use Getopt::Long;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $read;
my $write;
my $file;
my $help;

&GetOptions(
  'dbhost=s' => \$dbhost,
  'dbname=s' => \$dbname,
  'dbuser=s' => \$dbuser,
  'dbpass=s' => \$dbpass,
  'dbport=s' => \$dbport,
  'file=s'   => \$file,
  'read!'    => \$read,
  'write!'   => \$write,
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

my $ea = $db->get_EventAdaptor;
my $ra = $db->get_WorkflowAdaptor;

throw(
"Can't setup a workflow table or write a workflow config file without a file specified"
) unless ($file);

if ($read) {
  my $workflows = parse_file( $file, $ea );
  foreach my $workflow (@$workflows) {
    $ra->store($workflow);
  }
}
elsif ($write) {
  if ( -e $file ) {
    print STDERR $file . " exists this process will over write it\n";
  }
  open( my $fh, ">", $file ) or throw( "Failed to open " . $file );
  my $workflows = $ra->fetch_all;
  foreach my $workflow (@$workflows) {
    print $fh "[" . $workflow->goal_event->name . "]\n";
    foreach my $condition ( @{ $workflow->conditions } ) {
      print $fh "condition=" . $condition->name . "\n";
    }
  }
  close($fh);
}
else {
  throw(
    "Don't know what to run in workflow_setup neither read nor write has be "
      . "specified" );
}

sub read_file {
  my ($file) = @_;

  open( my $fh, "<", $file ) or throw "Couldn't open file $file";
  my %config;
  my $header;
  while (<$fh>) {
    chomp();
    next if ( /^\s$/ || /^\#/ );
    if (/^\[(.*)\]\s*$/) {    # $1 will be the header name, without the []
      $header = $1;
      if ( $config{$header} ) {
        warning( $header
            . " already exists in the file, you are now adding extra "
            . " conditions to it" );
      }
      $config{$header} = [] unless ( $config{$header} );
    }
    if (/^([^=\s]+)\s*=\s*(.+?)\s*$/) {    #key=value this defined pair
      my $key =
        lc($1);    # keys stored as all lowercase, values have case preserved
      my $value = $2;
      if ( length($header) == 0 ) {
        throw("Found key/value pair $key/$value outside stanza");
      }
      if ( $key eq 'condition' ) {
        push( @{ $config{$header} }, $value );
      }
      else {
        warning(
          "Don't recognise " . $key . " so doing nothing with " . $value );
      }
    }
  }
  close($fh);

  return %config;
}

sub parse_file {
  my ( $file, $ea ) = @_;

  my %config = read_file($file);

  my @workflows;
  foreach my $goal ( keys(%config) ) {
    my $workflow;
    eval {
      my $event = $ea->fetch_by_name($goal);
      
      throw( "Failed to find event for " . $goal . " in " . $ea->dbc->dbname )
        unless ($event);
      
      $workflow = ReseqTrack::Workflow->new( -goal_event => $event, );
      
      my $conditions = $config{$goal};
      
      foreach my $condition_name (@$conditions) {
        my $condition_event = $ea->fetch_by_name($condition_name);
        throw("Failed to find event for "
            . $condition_name . " in "
            . $ea->dbc->dbname )
          unless ($condition_event);
        $workflow->conditions($condition_event);
      }
    };
    if ($@) {
      throw("Failed to create a Workflow object for $goal $@");
    }
    push( @workflows, $workflow );
  }
  return \@workflows;
}

sub perldocs {
  exec( 'perldoc', $0 );
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/workflow/load_workflow_from_conf.pl

=head1 SYNOPSIS

This script can load workflow objects into a ReseqTrack database based on a config 
file in the format

[goal_event_name]

condition=conditional_event_name

A single goal can have multiple conditions if required.

=head1 OPTIONS

-dbhost, host name for database instance

-dbname, name of database on instance

-dbuser, name of user

-dbpass, password string if required

-dbport, the port number of the database instance

-file, name of file to read in or write out

-read, binary flag to indicate the file specified by -file should be parsed and loaded into the database

-write, binary flag to indicate the file specified should be written in the given format based on the workflow objects already in the database

-help, binary flag to indicate the help should be printed

=head1 Examples

perl ReseqTrack/scripts/workflow/load_workflow_from_conf.pl -dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database -file /path/to/workflow.conf -read

This will load the workflow objects specifed in workflow.conf

perl ReseqTrack/scripts/workflow/load_workflow_from_conf.pl -dbhost mysql-host -dbuser ro_user -dbport 4197 -dbname my_database -file /path/to/workflow.conf -write

This dumps the existing workflow objects into a file of the given name

The write process will overwrite any existing file with the given name

=head2 Other useful scripts

ReseqTrack/event/load_event_from_conf.pl is another script which works in a very similar manner but for event objects

=cut

