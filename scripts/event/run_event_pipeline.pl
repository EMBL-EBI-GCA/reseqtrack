#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Intersection;
use ReseqTrack::Tools::PipelineUtils;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Job;
use Getopt::Long;
use List::Util qw (shuffle);
use File::Basename;

$| = 1;

my $term_sig =  0;

$SIG{TERM} = \&termhandler;
$SIG{INT} = \&termhandler;


my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my @event_names;
my @allowed_inputs;
my $retry_max = 4;
my $once = 0;
my $submission_limit;
my $submission_total = 25;
my $num_output_dirs = 100;
my $runner;
my $max_number_of_jobs = 700000;
my $sleep = 360;
my $help;
my $batch_submission_path = 'ReseqTrack::Tools::BatchSubmission';
my $module_name = 'LSF';
my $verbose;
my $test;
my $override_lock;
my $submit_all_jobs = 0;
my $inputs_file;

&GetOptions( 
	    'dbhost=s'      => \$dbhost,
	    'dbname=s'      => \$dbname,
	    'dbuser=s'      => \$dbuser,
	    'dbpass=s'      => \$dbpass,
	    'dbport=s'      => \$dbport,
	    'name=s@'       => \@event_names,
	    'inputs=s@'     => \@allowed_inputs,
	    'max_retry=s'   => \$retry_max,
	    'submission_limit!' => \$submission_limit,
	    'submission_total=s' => \$submission_total,
	    'once!'         => \$once,
	    'output_dir_num=s' => \$num_output_dirs,
	    'runner=s'      => \$runner,
	    'max_numbers_of_jobs=s' => \$max_number_of_jobs,
	    'help!'         => \$help,
	    'batch_submission_path=s' => \$batch_submission_path,
	    'module_name:s' => \$module_name,
	    'sleep=s'       => \$sleep,
	    'verbose!'      => \$verbose,
	    'test!'         =>\$test,
	    'override_lock!'=>\$override_lock,
	    'submit_all_jobs!'=>\$submit_all_jobs,
	    'inputs_file=s' => \$inputs_file,
	   );


if ($help) {
    perldocs();
}
if ($test) {
    $once = 1;
}

if (!$runner){
	# guess the runner script using the default name and the dir of the running script
	my ($script_name,$script_path) = fileparse($0);
	$runner = $script_path."runner.pl";
	print "runner was not specified, will try to use $runner $/ " if ($verbose);
}

throw($runner." runner script must exist") unless(-e $runner);

print "Starting to run pipeline based on ".$dbname."\n" if($verbose);
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
					   -host => $dbhost,
					   -user => $dbuser,
					   -port => $dbport,
					   -dbname => $dbname,
					   -pass => $dbpass,
					  );
my $ja = $db->get_JobAdaptor;
my $meta_adaptor = $db->get_MetaAdaptor;

my $lock_is_set = 0;
eval{
    is_locked("pipeline.lock", $meta_adaptor);
};
if ($@) {
    if (! $override_lock) {
        throw("Can't run, there is a problem with the locking system \n");
    }
    else {
        print "System is locked, but override_lock flag is set\n";
    }
}
else {
    print "Locking the pipeline\n";
    create_lock_string("pipeline.lock", $meta_adaptor);
    $lock_is_set = 1;
}




if ($inputs_file) {
    open my $FILE, '<', $inputs_file
        or throw "can't open $inputs_file: $!";
    LINE:
    while (my $line = <$FILE>) {
        push(@allowed_inputs, $line);
    }
    close $FILE;
}

my %allowed_inputs;
INPUT:
foreach my $allowed_input (@allowed_inputs) {
    my ($table_name, $type, @inputs) = split(/[:\s]+/, $allowed_input, 3);
    next INPUT if (!$table_name);
    if (!$allowed_inputs{$table_name}) {
        $allowed_inputs{$table_name} = {};
    }
    next INPUT if (!$type);
    if (!$allowed_inputs{$table_name}{$type}) {
        $allowed_inputs{$table_name}{$type} = {};
    }
    foreach my $input (@inputs) {
        $allowed_inputs{$table_name}{$type}{$input} = 1;
    }
	if (lc($table_name) eq 'input_string') {
		my $input = join ' ', @inputs;
		$allowed_inputs{$table_name}{$type}{$input} = 1;
	}
}

my %event_name_hash;
foreach my $logic (@event_names) {
    $event_name_hash{$logic} = 1;
}

my %batch_sub_args;
$batch_sub_args{'-max_job_number'} = $max_number_of_jobs;
$batch_sub_args{'-sleep'} = $sleep;

my $batch_submission_module = $batch_submission_path."::".$module_name;
my $batch_submission_object = setup_batch_submission_system($batch_submission_module,
                                                            \%batch_sub_args);

print "Using ".$batch_submission_module." to submit jobs\n" if($verbose);
my ($events, $type_event_hash) = setup_event_objects($db, undef, \%event_name_hash);
print "Have ".@$events." events to consider\n" if($verbose);
my $workflow_hash = setup_workflow($db);
print "Have ".keys(%$workflow_hash)." workflow objects to use\n" if($verbose);
my $total_submission_count = 0;

my $loop_count = 1;
LOOP:
while(1){
    print "\n*****Starting loop $loop_count*****\n" if($verbose);
    $loop_count++;
    my $submitted_count = 0;
    my %submitted_hash;
    my %input_hash;
    EVENT:
    foreach my $event (@$events) {
        print "\nConsidering event ".$event->name."\n" if($verbose);
        my $table_name = $event->table_name;
        my $type = $event->type;

        my $check_allowed_inputs = 0;
        if (keys %allowed_inputs) {
            if (! exists $allowed_inputs{$table_name}) {
                print "Inputs not allowed for table ".$table_name."\n" if ($verbose);
                next EVENT;
            }
            if (keys %{$allowed_inputs{$table_name}}) {
                if (! exists $allowed_inputs{$table_name}{$type}) {
                    print "Inputs not allowed for table ".$table_name
                            ." and type ".$type."\n" if ($verbose);
                    next EVENT;
                }
                if (keys %{$allowed_inputs{$table_name}{$type}}) {
                    $check_allowed_inputs = 1;
                }
            }
        }

        my $existing_jobs = $ja->fetch_by_event($event);
        my $inputs = get_incomplete_inputs($db, $event);
        print "There are ".@$inputs." inputs that are not complete\n" if ($verbose);
        next EVENT if (!@$inputs);

        $inputs = check_existing_jobs($inputs, $existing_jobs, $retry_max, $verbose);
        print "There are ".@$inputs." inputs after filtering for existing jobs\n" if ($verbose);
        next EVENT if (!@$inputs);

        if ($check_allowed_inputs) {
            my @filtered_inputs = grep {$allowed_inputs{$table_name}{$type}{$_}} @$inputs;
            $inputs = \@filtered_inputs;
            print "There are ".@$inputs." allowed inputs\n" if ($verbose);
            next EVENT if (!@$inputs);
        }

        my $workflow = $workflow_hash->{$event->name};
        if ($workflow) {
            $inputs = check_workflow($workflow, $inputs, \%input_hash, $db, $verbose);
            print "There are ".@$inputs." inputs after checking the workflow\n" if ($verbose);
            next EVENT if (!@$inputs);
        }

        if ($term_sig) {
            print "Have term_sig so exiting loop\n";
            last LOOP;
        }

        next EVENT if (!@$inputs);
        print "Have ".@$inputs." inputs to submit for ".$event->name."\n";

        if ($test) {
            my $submission_cmd =
                    create_test_submission_cmd($db, $batch_submission_object, $runner, $event, scalar @$inputs);
            print $submission_cmd, "\n";
            foreach my $input (@$inputs) {
                my $cmd = create_event_commandline($event, $input);
                print $cmd, "\n";
            }
            next EVENT;
        }
    
        my $num_lsf_jobs = scalar grep {/RUNNING/ || /SUBMITTED/ || /WAITING/} map {$_->current_status} @$existing_jobs;
        my @jobs_to_submit;
        INPUT:
        foreach my $input (shuffle @$inputs) {
            last INPUT if ($submission_limit && $total_submission_count + @jobs_to_submit >= $submission_total);
            if (!$submit_all_jobs && $num_lsf_jobs) {
              last INPUT if (@jobs_to_submit + $num_lsf_jobs >= $event->job_slot_limit);
            }
            my $job = create_job_object($event, $input, $db, $num_output_dirs, $retry_max);
            push(@jobs_to_submit, $job) if($job);
        }
  
        SUBMIT:
        while (scalar @jobs_to_submit) {
            my $job_slots_available = $batch_submission_object->job_slots_available();
    
            my @jobs;
            JOB:
            while (scalar @jobs_to_submit) {
                last JOB if (scalar @jobs >= $job_slots_available);
    
                my $job = shift @jobs_to_submit;
                push(@jobs, $job);
                $total_submission_count ++;
                $submitted_count ++;
            }

            my $cmds = create_submission_cmds(\@jobs, $db, $batch_submission_object, $runner, $event);
            foreach my $cmd (keys %$cmds) {
                foreach my $job (@{$cmds->{$cmd}}) {
                    $job->current_status("SUBMITTED");
                    $ja->set_status($job);
                    print "Submitting ".$job->dbID."\n" if($verbose);
                    if ($submitted_hash{$job->dbID}) {
                        throw("Have already submitted ".$job->dbID);
                    } else {
                        $submitted_hash{$job->dbID} = 1;
                    }
                }
                eval{
                    system($cmd);
                    print "CMD ".$cmd."\n" if($verbose);
                };
                if ($@) {
                    print STDERR "Failed to run ".$cmd." $@\n";
                    foreach my $job (@{$cmds->{$cmd}}) {
                        $job->current_status("FAILED");
                        $ja->set_status($job);
                    }
                }
            }

        } # end of SUBMIT loop
    } # end of EVENT loop


    #check if pipeline is finished
    if ($term_sig) {
        print "Have term_sig so exiting loop\n";
        last LOOP;
    }
    if ($once) {
        print "\'once\' flag is set so exiting loop\n";
        last LOOP;
    }
    if ($submission_limit && $total_submission_count >= $submission_total) {
        print "reached submission limit so exiting loop\n";
        last LOOP;
    }
    if (!$submitted_count && check_if_done($db, $retry_max)) {
        print "pipeline complete so exiting loop\n";
        last LOOP;
    }

} # end of LOOP


if ($lock_is_set) {
    print "unlocking the pipeline\n";
    delete_lock_string("pipeline.lock", $meta_adaptor);
}
  
sub termhandler {
    $term_sig = 1;
}

sub perldocs{
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

-inputs, this specifies which input strings will be considered.
Acceptable values are table_name (e.g. -inputs=collection),
table_name:input_type (e.g. -inputs=collection:FILTERED_FASTQ)
or table_name:input_type:input (e.g. -inputs=collection:FILTERED_FASTQ:ERR000001)
This can appear on the command line multiple times.  All input strings will be considered if none are specified.

-inputs_file, the path to a file containing inputs to be considered.
For each line in the file, first field is table_name, second field is input_type (optional), third field is input_string (optional).
Fields are separated by whitespace.  Lines starting with # are ignored.  All input strings will be considered if none are specified.

-max_retry, this is the maximum number of times a single job will be retried if it fails, this defaults to 4

-submission_limit, this indicates the script should only loop once and should limit
the number of jobs submitted, by default this is 25

-submission_total, this specified the number of jobs -submission_limit should be allowed to be submitted

-once, this flag means the script only runs through its loop once so no failed jobs will be retried nor no dependancies which aren't already satified reached

-test, this flag means jobs are not submitted to the batch submission system.  Instead, the script will print the program command line.

-submit_all_jobs, this flag is relevant when a large array of jobs has already been submitted to lsf.  When this flag is set, any new jobs will be submitted even if this results in the number of running jobs exceeding the value of job_slot_limit.  Use with caution.

-override_lock, this flag means the script will continue to run even if the pipeline is locked

-output_dir_number, job output files and farm log files are split across in the event specified output_dir; by default this is 100.

-runner, this is a path to the runner script which is used to run submitted jobs. Generally ReseqTrack/scripts/event/runner.pl

-max_number_of_jobs, the EBI farm does not allow more than 700000 pending jobs in the system at once, once this limit is reached the script sleeps for 5 minutes and checks again, continuing this until there are less than 700000 pending jobs, if you wish to change this limit you need to alter this argument

-batch_submission_path, root perl path for batch submission module; this defaults to
ReseqTrack::Tools::BatchSubmission

-module_name, name of batch submission module; this defaults to LSF

-verbose, make progress print statements

-help, binary flag to indicate the help should be printed

=head1 Examples

perl ReseqTrack/scripts/event/run_event_pipeline.pl -dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database -runner ReseqTrack/scripts/event/runner.pl -name my_event -once

This would only submit jobs for the event, my_event and would only loop once through the data.

perl ReseqTrack/scripts/event/run_event_pipeline.pl -dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database -runner ReseqTrack/scripts/event/runner.pl 

This would submit jobs for all incomplete events.  It will continue to run through the loop until all running jobs have completed and no more jobs can be submitted.


=head1 Other useful scripts

ReseqTrack/scripts/event/montior_event_pipeline.pl may also be useful as it summarizes what input is available, what is running and what has finished.

=cut

