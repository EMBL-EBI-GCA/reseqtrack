#!/sw/arch/bin/perl -w

use strict;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Intersection;
use ReseqTrack::Tools::PipelineUtils;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Job;
use Getopt::Long;
use File::Copy;
use File::Path;
use File::Basename;
use Benchmark;
$| = 1;

my $term_sig = 0;

$SIG{TERM} = \&termhandler;
$SIG{INT}  = \&termhandler;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my @event_names;
my @input_types;
my $retry_max = 4;
my $once      = 0;
my $submission_limit;
my $submission_total = 25;
my $num_output_dirs  = 10;
my $runner;
my $max_number_of_jobs = 10000;
my $sleep              = 360;
my $help;
my $batch_submission_path = 'ReseqTrack::Tools::BatchSubmission';
my $module_name           = 'LSF';
my $verbose;

&GetOptions(
    'dbhost=s'                => \$dbhost,
    'dbname=s'                => \$dbname,
    'dbuser=s'                => \$dbuser,
    'dbpass=s'                => \$dbpass,
    'dbport=s'                => \$dbport,
    'name=s@'                 => \@event_names,
    'type=s@'                 => \@input_types,
    'max_retry=s'             => \$retry_max,
    'submission_limit!'       => \$submission_limit,
    'submission_total=s'      => \$submission_total,
    'once!'                   => \$once,
    'output_dir_num=s'        => \$num_output_dirs,
    'runner=s'                => \$runner,
    'max_numbers_of_jobs=s'   => \$max_number_of_jobs,
    'help!'                   => \$help,
    'batch_submission_path=s' => \$batch_submission_path,
    'module_name:s'           => \$module_name,
    'sleep=s'                 => \$sleep,
    'verbose!'                => \$verbose,
);

if ($help) {
    perldocs();
}

print "Starting to run pipeline based on " . $dbname . "\n" if ($verbose);
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);

my $meta_adaptor = $db->get_MetaAdaptor;
eval { is_locked("pipeline.lock", $meta_adaptor); };
if ($@) {
    throw("Can't run, there is a problem with the locking system \n" . "$@\n");
}
create_lock_string("pipeline.lock", $meta_adaptor);
my $ja = $db->get_JobAdaptor;

throw($runner . " runner script must exist") unless (-e $runner);
my %allowed_input_types;
foreach my $type (@input_types) {
    $allowed_input_types{$type} = 1;
}
my %event_name_hash;
foreach my $logic (@event_names) {
    $event_name_hash{$logic} = 1;
}

my %batch_sub_args;
$batch_sub_args{'-max_job_number'} = $max_number_of_jobs;
$batch_sub_args{'-sleep'}          = $sleep;

my $batch_submission_module = $batch_submission_path . "::" . $module_name;
my $batch_submission_object =
  setup_batch_submission_system($batch_submission_module, \%batch_sub_args);

print "Using " . $batch_submission_module . " to submit jobs\n" if ($verbose);
my ($events, $type_event_hash) =
  setup_event_objects($db, undef, \%event_name_hash);
print "Have " . @$events . " events to consider\n" if ($verbose);
my $workflow_hash = setup_workflow($db);
print "Have " . keys(%$workflow_hash) . " workflow objects to use\n"
  if ($verbose);
my $total_submission_count = 0;
my $done                   = 0;

my %event_hash;
foreach my $event (@$events) {
    $event_hash{$event->name} = $event;
}
my $loop_count = 1;
LOOP: while (1) {
    my $submitted_count = 0;
    print "\n*****Starting loop $loop_count*****\n" if ($verbose);
    $loop_count++;
    my %submitted_hash;
    my $input_hash = get_inputs($db, $events);
    my %event_to_submit;
  TABLE_NAME: foreach my $table_name (keys(%$input_hash)) {
        print "TABLE NAME " . $table_name . "\n" if ($verbose);
        my %type_hash = %{$input_hash->{$table_name}};
      INPUT_TYPE: foreach my $input_type (keys(%type_hash)) {
            if (keys(%allowed_input_types)) {
                next INPUT_TYPE unless ($allowed_input_types{$input_type});
            }
            print "INPUT TYPE " . $input_type . "\n" if ($verbose);
            my $inputs = $type_hash{$input_type};
            my $events = $type_event_hash->{$input_type};
            foreach my $input (@$inputs) {
                if ($term_sig) {
                    print "Have sigterm\n";
                    $done = 1;
                    last INPUT_TYPE;
                }
              EVENT: foreach my $event (@$events) {
                    next EVENT if ($event->table_name ne $table_name);
                    my $workflow = $workflow_hash->{$event->name};
                    next
                      if (@event_names && !$event_name_hash{$event->name});
                    print "Checking if "
                      . $event->name
                      . " can run on "
                      . $input . "\n"
                      if ($verbose);
                    if (
                        can_run(
                            $input,      $event,     $db, $workflow,
                            $input_hash, $retry_max, $verbose
                        )
                      )
                    {
                        if (!$event_to_submit{$event->name}) {
                            $event_to_submit{$event->name} = [];
                        }
                        push(@{$event_to_submit{$event->name}}, $input);
                    } else {
                        print "Cant run\n" if ($verbose);
                    }
                }
            }
        }
        if ($done) {
            last TABLE_NAME;
        }
    }

  EVENT_NAME: foreach my $name (keys(%event_to_submit)) {
        print "\nConsidering event " . $name . "\n" if ($verbose);
        if ($term_sig) {
            $done = 1;
            last EVENT_NAME;
        }
        if ($done) {
            last EVENT_NAME;
        }
        my $inputs_to_submit = $event_to_submit{$name};
        my $event            = $event_hash{$name};
        my @jobs;
        my %input_strings;
        foreach my $input (@$inputs_to_submit) {
            my $job = create_job_object($event, $input, $db, $num_output_dirs,
                $retry_max);
            push(@jobs, $job) if ($job);
        }
        my $cmds =
          create_submission_cmds(\@jobs, $db, $batch_submission_object, $runner,
            $event);
      CMD: foreach my $cmd (keys(%$cmds)) {
            if ($submission_limit) {
                if ($total_submission_count >= $submission_total) {
                    $done = 1;
                    last CMD;
                }
            }
            my $return_value =
              has_to_many_jobs($max_number_of_jobs, $batch_submission_object,
                $sleep);
            my $jobs = $cmds->{$cmd};
            foreach my $job (@$jobs) {
                $job->current_status("SUBMITTED");
                $ja->set_status($job);
                print "Submitting " . $job->dbID . "\n" if ($verbose);
                if ($submitted_hash{$job->dbID}) {
                    throw("Have already submitted " . $job->dbID);
                } else {
                    $submitted_hash{$job->dbID} = 1;
                }
            }
            eval {
                system($cmd);
                print $cmd. "\n";
                $total_submission_count++;
                $submitted_count++;
            };
            if ($@) {
                print STDERR "Failed to run " . $cmd . " $@\n";
                foreach my $job (@$jobs) {
                    $job->current_status("FAILED");
                    $ja->set_status($job);
                }
                next CMD;
            }
        }
    }

    #check if pipeline is finished

    if ($once || $done) {
        last LOOP;
    } else {
        if (check_if_done($db, $retry_max)) {
            last LOOP unless ($submitted_count >= 1);
        }
    }
}

delete_lock_string("pipeline.lock", $meta_adaptor);

sub termhandler {
    $term_sig = 1;
}

sub perldocs {
    exec('perldoc', $0);
    exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/event/run_event_pipeline.pl

=head1 SYNOPSIS

This script uses information present in the event, event_complete, job and workflow
tables to deterime what events can be run with what input on the farm. 

=head1 OPTIONS

-dbhost, host name for database instance

-dbname, name of database on instance

-dbuser, name of user

-dbpass, password string if required

-dbport, the port number of the database instance

-name, this is the name of an entry in the event table, if specified only the 
event names given will be considered for jobs to be submitted, this can appear on the 
commandline multiple times

-type, this specifies the type of input which will be considered, once specified only types on the list will be used, this can also appear on the commandline multiple times

-max_retry, this is the maximum number of times a single job will be retried if it fails, this defaults to 4

-submission_limit, this indicates the script should only loop once and should limit
the number of jobs submitted, by default this is 25

-submission_total, this specified the number of jobs -submission_limit should be allowed to be submitted

-once, this means the script only runs through its loop once so no failed jobs will be retried nor no dependancies which aren't already satified reached

-run, this defaults to on and means system calls are actually made on the bsub commands

-output_dir_number, STDERR and STDOUT files are split across in the event specified output_dir, by default this is 10, 0-9 if you are going to run significantlly more than 10000 jobs

-runner, this is a path to the runner script which is used to run submitted jobs. Generally ReseqTrack/scripts/event/run

-max_number_of_jobs, the EBI farm does not allow more than 10000 pending jobs in the system at once, once this limit is reached the script sleeps for 5 minutes and checks again, continuing this untill there are less than 10000 pending jobs, if you wish to change this limit you need to alter this argument

-batch_submission_path, root perl path for batch submission module, this defaults to
ReseqTrack::Tools::BatchSubmission

-module_name, name of batch submission module, this defaults to LSF

-verbose, make progress print statements

-help, binary flag to indicate the help should be printed

=head1 Examples

perl ReseqTrack/scripts/event/run_event_pipeline.pl -dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database -runner ReseqTrack/scripts/event/runner.pl -name my_event -once

This would only submit jobs for the event, my_event and would only loop once through the data.

perl ReseqTrack/scripts/event/run_event_pipeline.pl -dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database -runner ReseqTrack/scripts/event/runner.pl 

This would deterime what events it could run on the basis of the information in the database. It would also keep running over the loop untill it recieves a sigterm 


=head1 Other useful scripts

ReseqTrack/scripts/event/montior_event_pipeline.pl may also be useful as it summarizes what input is available, what is running and what has finished.

=cut
