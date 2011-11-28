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
use Data::Dumper;
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
my $retry_max = 1;
my $once = 0;
my $submission_limit;
my $submission_total = 25;
my $num_output_dirs = 100;
my $runner;
my $max_number_of_jobs = 10000;
my $sleep = 360;
my $help;
my $batch_submission_path = 'ReseqTrack::Tools::BatchSubmission';
my $module_name = 'LSF';
my $verbose;
my $submit_hash;
my $test;
my $override_lock;
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
	    'inputs_file=s' => \$inputs_file,
	   );




if ($help) {
    perldocs();
}
if ($test) {
    $once = 1;
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
    my ($table_name, $type, @inputs) = split(/[:\s]+/, $allowed_input);
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

my %event_hash;
foreach my $event (@$events) {
    $event_hash{$event->name} = $event;
}

my $loop_count = 1;
LOOP:
while(1){
    print "\n*****Starting loop $loop_count*****\n" if($verbose);
    $loop_count++;
    my $submitted_count = 0;
    my %submitted_hash;
    my $input_hash = get_inputs($db, $events);
    my %event_to_submit;
    TABLE_NAME:
    foreach my $table_name(keys(%$input_hash)){
        next TABLE_NAME if (keys %allowed_inputs && !$allowed_inputs{$table_name});
        print "TABLE NAME ".$table_name."\n" if($verbose);
        my %type_hash = %{$input_hash->{$table_name}};
        INPUT_TYPE:
        foreach my $input_type(keys(%type_hash)){
            next INPUT_TYPE if ($allowed_inputs{$table_name} && keys %{$allowed_inputs{$table_name}}
                                && !$allowed_inputs{$table_name}{$input_type});
            my $inputs = $type_hash{$input_type};
            print "INPUT TYPE ".$input_type."\n" if($verbose);
            my $events_for_input_type = $type_event_hash->{$input_type};
            INPUT:
            foreach my $input (@$inputs) {
                next INPUT if ($allowed_inputs{$table_name}{$input_type} &&
                                keys %{$allowed_inputs{$table_name}{$input_type}}
                                && !$allowed_inputs{$table_name}{$input_type}{$input});
                EVENT:
                foreach my $event(@$events_for_input_type){
                    next EVENT if($event->table_name ne $table_name);
                    next EVENT if(@event_names && !$event_name_hash{$event->name});
                    my $workflow = $workflow_hash->{$event->name};
                    print "Checking if ".$event->name." can run on ".$input."\n" if($verbose);
                    if (can_run($input, $event, $db, $workflow, $input_hash, $retry_max, $verbose)) {
                        if (!$event_to_submit{$event->name}) {
                            $event_to_submit{$event->name} = [];
                        }
                        push(@{$event_to_submit{$event->name}}, $input);
                    } else {
                        print "Cant run\n" if($verbose);
                    }


                } # end of EVENT loop
            } # end of INPUT loop
        } # end of INPUT_TYPE loop
    } # end of TABLE_NAME loop

    EVENT_NAME:
    foreach my $name(keys(%event_to_submit)){
        print "\nConsidering event ".$name."\n";# if($verbose);

        my $inputs_to_submit = $event_to_submit{$name};
        my $event = $event_hash{$name};
        my @jobs_to_submit;
        my %input_strings;
    
        print "Have ".@$inputs_to_submit." inputs to submit\n"; # if($verbose);
    
        if ($test) {
            foreach my $input (@$inputs_to_submit) {
                my $cmd = create_event_commandline($event, $input);
                print $cmd, "\n";
            }
            next EVENT_NAME;
        }
    
        foreach my $input (@$inputs_to_submit) {
            my $job = create_job_object($event, $input, $db, $num_output_dirs, $retry_max);
            push(@jobs_to_submit, $job) if($job);
        }
  
        SUBMIT:
        while (scalar @jobs_to_submit) {
            if ($submission_limit && $total_submission_count >= $submission_total) {
                print "reached submission limit so exiting loop\n";
                last LOOP;
            }
            my $job_slots_available = $batch_submission_object->job_slots_available();
    
            my @jobs;
            JOB:
            while (scalar @jobs_to_submit) {
                last JOB if (scalar @jobs >= $job_slots_available);
                last JOB if ($submission_limit && $total_submission_count >= $submission_total);
    
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
    } # end of EVENT_NAME loop


    #check if pipeline is finished
    if ($term_sig) {
        print "Have term_sig so exiting loop\n";
        last LOOP;
    }
    if ($once) {
        print "\'once\' flag is set so exiting loop\n";
        last LOOP;
    }
    if (check_if_done($db, $retry_max) && !$submitted_count) {
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

-override_lock, this flag means the script will continue to run even if the pipeline is locked

-output_dir_number, job output files and farm log files are split across in the event specified output_dir; by default this is 100.

-runner, this is a path to the runner script which is used to run submitted jobs. Generally ReseqTrack/scripts/event/runner.pl

-max_number_of_jobs, the EBI farm does not allow more than 10000 pending jobs in the system at once, once this limit is reached the script sleeps for 5 minutes and checks again, continuing this untill there are less than 10000 pending jobs, if you wish to change this limit you need to alter this argument

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

