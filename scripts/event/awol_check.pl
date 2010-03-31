#!/sw/arch/bin/perl5.8.7 -w

use strict;
use Getopt::Long;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $user;
my $help;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'user=s' => \$user,
  'help!' => \$help,
    );

if ($help) {
  perldocs();
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );

$user = $ENV{USER} if(!$user);
my $ja = $db->get_JobAdaptor;
my $jobs = $ja->fetch_all;


my $hash = run_bjobs($user);
my $count = 0;
foreach my $job(@$jobs){
  next if($job->current_status =~ /^FAIL/);
  next if($job->current_status =~ /AWOL/);
  next if($job->current_status eq 'CREATED');
  next if(!$job->submission_id || $job->submission_id == 0);
  unless($hash->{$job->submission_id}){
    $job->current_status("AWOL");
    $ja->set_status($job);
  }else{
    #print $job->dbID." ".$job->submission_id." ".$hash->{$job->submission_id}."\n";
  }
}

sub run_bjobs{
  my ($user) = @_;
  my $cmd = "bjobs -w -u $user";
  print $cmd."\n";
  open(CMD, $cmd." |") or throw("Failed to run ".$cmd);
  my %hash;
  while(<CMD>){
    #print;
    chomp;
    my @values = split;
    next unless($values[0] =~ /\d+/);
    my $submission_id = $values[0];
    my $job_name = $values[6];
    $hash{$submission_id} = 1;
  }
  return \%hash;
}

sub perldoc{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/files/load_event_from_conf.pl

=head1 SYNOPSIS

=head2 OPTIONS

   -dbhost, the name of the mysql-host

   -dbname, the name of the mysql database

   -dbuser, the name of the mysql user

   -dbpass, the database password if appropriate

   -dbport, the port the mysql instance is running on, this defaults to 4197 the 
   standard port for mysql-g1kdcc.ebi.ac.uk

   -user, this allows the script runner to specify a different lsf user to themselves
    without this defined it will only get jobs which are owned by the script runner

   -help, binary flag to print out perldocs

=head1 Example:

perl ReseqTrack/scripts/event/awol_check.pl -dbhost host -dbuser rwuser -dbpass **** -dbport 4197 -dbname my_database

=head2 Other useful scripts

ReseqTrack/scripts/event/monitor_event_pipeline.pl

=cut
