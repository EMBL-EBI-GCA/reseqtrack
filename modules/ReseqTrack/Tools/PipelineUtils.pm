=pod

=head1 NAME

ReseqTrack::Tools::PipelineUtils;

=head1 SYNOPSIS

This is a collection of methods useful for running the Event Pipeline section of
the ReseqTrack schema

=head1 Example

use ReseqTrack::Tools::PipelineUtils qw (get_inputs);

my $input_hash = get_inputs($db, $events);

=cut

package ReseqTrack::Tools::PipelineUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::Job;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use File::Path;
use File::Basename;
use Sys::Hostname;
use Socket;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);

@EXPORT = qw(get_inputs get_all_inputs create_event_commandline
          get_random_input_string setup_event_objects setup_batch_submission_system
          can_run check_workflow create_job_object create_submission_cmds
          has_to_many_jobs setup_workflow check_if_done create_job_output_stem);

=head2 get_all_inputs

  Arg [1]   : ReseqTrack::DBSQL::DBAdaptor
  Function  : create a list of all input strings
  Returntype: arrayref list of input strings
  Exceptions: n/a
  Example   : my $inputs = get_all_inputs($db);

=cut


sub get_all_inputs{
  my ($db) = @_;
  my $events = $db->get_EventAdaptor->fetch_all();
  my $input_hash = get_inputs($db, $events);
  my @inputs;
  foreach my $table_name(keys(%$input_hash)){
    foreach my $type(keys(%{$input_hash->{$table_name}})){
      my $strings = $input_hash->{$table_name}->{$type};
      push(@inputs, @$strings);
    }
  }
  return \@inputs;
}



=head2 get_inputs

  Arg [1]   : ReseqTrack::DBSQL::DBAdaptor;
  Arg [2]   : arrayref of ReseqTrack::Event objects
  Function  : If not passed events, it fetches them all from the database. The
  is uses the table names and types associated with each event to fetch input strings
  that each event should be run on
  Returntype: hashref, returns a hashref, this is a hash of hashes keyed first on 
  table name then on input type the values are arrayrefs of input strings
  Exceptions: throws if doesn't recognise a table name and how to fetch inputs from
  that table
  Example   : my $input_hash = get_inputs($db, $event); 

=cut



sub get_inputs{
  my ($db, $events) = @_;
  if(!$events){
    $events = $db->get_EventAdaptor->fetch_all();
  }

  my %tables_to_input;

  foreach my $event(@$events){
    $tables_to_input{$event->table_name}{$event->type} = [] 
        if(!$tables_to_input{$event->table_name}{$event->type});

    next if(@{$tables_to_input{$event->table_name}{$event->type}});

    if($event->table_name eq 'file'){
      my $inputs = get_file_inputs($db, $event->type);
      $tables_to_input{$event->table_name}{$event->type} = $inputs;
    }
    elsif($event->table_name eq 'collection'){
      my $inputs = get_collection_inputs($db, $event->type);
      $tables_to_input{$event->table_name}{$event->type} = $inputs;
    }
    elsif($event->table_name eq 'run_meta_info'){
      my $inputs = get_run_meta_info_inputs($db, $event->type);
      $tables_to_input{$event->table_name}{$event->type} = $inputs;
    }
    elsif($event->table_name eq 'input_string'){
      my $inputs = get_input_string_inputs($db, $event->type);
      $tables_to_input{$event->table_name}{$event->type} = $inputs;
    }
    else{
      throw("Don't know what sort of inputs to fetch for ".$event->name." ".
            $event->table_name);
    }
  }
  
  return \%tables_to_input;
}


=head2 get_XXXX_inputs

  Arg [1]   : ReseqTrack::DBSQL::DBAdaptor
  Arg [2]   : string, type
  Function  : to fetch input strings from given table based on given type
  Returntype: arrayref
  Exceptions:
  Example   : my $file_inputs = get_file_inputs($db, $type); 

=cut



sub get_file_inputs{
  my ($db, $type) = @_;
  my $fa = $db->get_FileAdaptor;
  my $files = $fa->fetch_by_type($type);
  my @inputs;
  foreach my $file(@$files){
    push(@inputs, $file->full_path);
  }
  return \@inputs;
}

sub get_collection_inputs{
  my ($db, $type) = @_;
  my $ca = $db->get_CollectionAdaptor;
  my $collections = $ca->fetch_by_type($type);
  my @inputs;
  foreach my $collection(@$collections){
    push(@inputs, $collection->name);
  }
  return \@inputs;
}

sub get_run_meta_info_inputs{
  my ($db, $type) = @_;
  my $ca = $db->get_RunMetaInfoAdaptor;
  my $run_meta_infos = $ca->fetch_all($type);
  my @inputs;
  foreach my $run_meta_info(@$run_meta_infos){
    push(@inputs, $run_meta_info->run_id);
  }
  return \@inputs;
}

sub get_input_string_inputs{
  my ($db, $type) = @_;
  my $ca = $db->get_InputStringAdaptor;
  my $run_meta_infos = $ca->fetch_by_type($type);
  my @inputs;
  foreach my $run_meta_info(@$run_meta_infos){
    push(@inputs, $run_meta_info->name);
  }
  return \@inputs;
}


=head2 create_event_commandline

  Arg [1]   : ReseqTrack::Event object
  Arg [2]   : input_string
  Function  : Create a command line from event table fields
  Returntype: string
  Example   : create_event_commandline($tevent, $file);

=cut

sub create_event_commandline{
    my ($event, $input_string) = @_;

    my $cmd = $event->program . " ";
    $cmd .= $event->options . " "    if ($event->options);
    $cmd .= $event->input_flag . " " if ($event->input_flag);
    $cmd .= $input_string;

    return $cmd;
}



=head2  get_random_input_string


  Arg [1]   : ReseqTrack::DBSQL::DBAdaptor
  Arg [2]   : ReseqTrack::Event
  Function  : generate a random input string for an event
  Returntype: string
  Exceptions: Must pass in an dbadaptor and an event object
  Example   : get_random_input_string ($db, $event);

=cut

sub get_random_input_string{
  my($db, $event) = @_;
  throw("Must pass get_randon_input_string a dbadaptor") 
      unless($db && $db->isa("ReseqTrack::DBSQL::DBAdaptor"));
  throw("Must pass get_random_input_string an event object") 
      unless($event && $event->isa("ReseqTrack::Event"));
  my $input_hash = get_inputs($db, [$event]);
  my $inputs = $input_hash->{$event->table_name}{$event->type};
  my $random_number = int ( rand ($#$inputs+1) );
  return $inputs->[$random_number];
}

=head2 setup_batch_submission_system

  Arg [1]   : string, path of batch submission module
  Function  : to include the batch submission module in the name space and
  to create the object for further use
  Returntype: batch submission object
  Exceptions: throws if it fails to require the given path
  Example   : 

=cut

sub setup_batch_submission_system{
  my ($batch_submission_module, $args) = @_;
  my $file = "$batch_submission_module.pm";
  $file =~ s{::}{/}g;
  eval {
    require "$file";
  };
  if($@){
    throw("ReseqTrack::Tools::EventPipeline::setup_batch_submission_system ".
          "Can't find $file [$@]");
  }
  my %constructor_args = %$args;
  my $object = $batch_submission_module->new(%constructor_args);
  return $object
}


=head2 setup_event_objects

  Arg [1]   : ReseqTrack::DBSQL::DBAdaptor
  Arg [2]   : arrayref of ReseqTrack::Event objects, if not supplied the list will
   be fetched from the database with fetch_all
  Arg [3]   : hashref, hash keyed on event names listing the allowed events. If 
  this is empty it is assumed all events are allowed
  Function  : process the events into more useful format for the running of the 
  pipeline
  Returntype: arrayref of ReseqTrack::Events and hashref keyed on event type pointing
  to an arrayref of events of that type 
  Exceptions: n/a
  Example   : my ($events, $event_type_hash) = setup_event_objects($db);

=cut


sub setup_event_objects{
  my ($db, $events, $allowed_event_names) = @_;
  $events = $db->get_EventAdaptor->fetch_all unless($events && @$events >= 1);
  my %event_hash;
  my @events;
  my %type_event_hash;
  foreach my $event(@$events){
    $event_hash{$event->name} = $event;
    $type_event_hash{$event->type} = [] if(!$type_event_hash{$event->type});
    push(@{$type_event_hash{$event->type}}, $event);
    if(keys(%$allowed_event_names)){
      push(@events, $event) if($allowed_event_names->{$event->name});
    }else{
      push(@events, $event);
    }
  }
  throw("ReseqTrack::Tools::PipelineUtils::setup_event_objects has no event ".
        "objects something has gone wrong") unless(\@events >= 1);
  foreach my $event(@events){
    throw("ReseqTrack::Tools::PipelineUtils::setup_event_objects ".$event->name.
          " doesn't have a table name defined") unless($event->table_name);
     throw("ReseqTrack::Tools::PipelineUtils::setup_event_objects ".$event->name.
          " doesn't have a type defined") unless($event->type);
  }
  return(\@events, \%type_event_hash);
}


=head2 setup_workflow

  Arg [1]   : ReseqTrack::DBSQL::DBAdaptor
  Arg [2]   : arrayref of ReseqTrack::Workflow objects
  Function  : turn array of workflow objects into a hash keyed on goal event name
  Returntype: hashref keyed on goal event name
  Exceptions: n/a
  Example   : 

=cut


sub setup_workflow{
  my ($db, $workflows) = @_;
  my $wfa = $db->get_WorkflowAdaptor;
  my %workflow_hash;
  $workflows = $wfa->fetch_all unless($workflows);
  foreach my $workflow(@$workflows){
    $workflow_hash{$workflow->goal_event->name} = $workflow;
  }
  return \%workflow_hash;
}



=head2 can_run

  Arg [1]   : string, input string
  Arg [2]   : ReseqTrack::Event
  Arg [3]   : ReseqTrack::DBSQL::DBAdaptor;
  Arg [4]   : ReseqTrack::Workflow
  Arg [5]   : hashref, hash of hashes, first keyed on table name, then event type
  values being arrays of input string
  Arg [6]   : int, maximum retry count for a job
  Function  : check if a particular event/input string combination can run
  Returntype: binary, 0/1 if true the job can run
  Exceptions: n/a
  Example   : if(can_run($input, $event, $db, $workflow, $input_hash, $max_retry));

=cut


sub can_run{
  my ($input, $event, $db, $workflow, $input_hash, $retry_max, $verbose) = @_;
  print "Trying to see if we can run ".$event->name." with ".$input."\n" 
      if($verbose);
  my $csa = $db->get_EventCompleteAdaptor;
  throw("Don't see to have an EventCompleteAdaptor") if(!$csa);
  
  my $ja = $db->get_JobAdaptor;
  
  if($csa->fetch_by_input_string_and_event($input, $event)){
    print $input." ".$event->name." is already complete\n" if($verbose);
    return 0;
  }else{
    my $workflow_conditions_fulfilled = 0;
    if($workflow){
      if(check_workflow($workflow, $event, $input, $input_hash, $csa, $verbose)){
        $workflow_conditions_fulfilled = 1;
      }
    }else{
      $workflow_conditions_fulfilled = 1;
    }
    unless($workflow_conditions_fulfilled){
      return 0;
    }
    my $job = $ja->fetch_by_name_and_input_string
        ($event->name, $input);
    print "No job in system can run ".$input." ".$event->name."\n" if($verbose && !$job);
    return 1 if(!$job);
    print "Have job ".$job->dbID." ".$job->current_status."\n" if($verbose);
    $job->current_status(undef, 1);
    if($job->current_status ne 'FAILED' &&
       $job->current_status ne 'AWOL' && 
       $job->current_status ne 'CREATED'){
      return 0;
    }elsif($job->current_status eq "FAILED" || 
           $job->current_status eq "AWOL"){
      return 0 if($job->retry_count && $job->retry_count >= $retry_max);
    }
    if($csa->fetch_by_input_string_and_event($input, $event)){
      print $input." ".$event->name." is already complete\n" if($verbose);
      return 0;
    }
    return 1;
  }
  throw("ReseqTrack::Tools::EventUtils:can_run, shouldn't of got this far ".
        "in can run");
}


=head2 check_workflow

  Arg [1]   : ReseqTrack::Workflow
  Arg [2]   : ReseqTrack::Event
  Arg [3]   : hashref. hash of hashes, first keyed on table name, then event type
  values being arrays of input string
  Arg [4]   : ReseqTrack::DBSQL::EventCompleteAdaptor
  Function  : check if conditions of workflow has been met
  Returntype: binary 0/1
  Exceptions: n/a
  Example   : if(check_workflow($workflow, $event, $input, $input_hash, $eca));

=cut



sub check_workflow{
  my ($workflow, $event, $input, $input_type_hash, $csa, $verbose) = @_;
  print "Checking workflow\n" if($verbose);
  my $conditions = $workflow->conditions;
  unless($conditions && $conditions >= 1){
    print "There are no conditions for workflow with ".
        $workflow->goal_event->name."\n" if($verbose);
    return 1;
  }else{
    print "There are ".@$conditions." associated with workflow for ".
        $workflow->goal_event->name."\n" if($verbose);
  }
  #Checking conditions
  print "Checking conditions for ".$event->name."\n" if($verbose);
  foreach my $condition(@$conditions){
    print "Looking at ".$condition->name."\n" if($verbose);
    my $checked = 0;
    if($condition->table_name eq $event->table_name){
      if($condition->type eq $event->type){
        print "table names and types match\n" if($verbose);
        $checked = 1;
        unless($csa->fetch_by_input_string_and_event($input, $condition)){
          print $condition->name." isn't complete on ".$input." can't run ".$event->name."\n"
              if($verbose);
          return 0;
        }
      }
    }
    unless($checked){
      my $completed_strings = $csa->fetch_by_event($condition);
      my $total_inputs = 
          $input_type_hash->{$condition->table_name}->{$condition->type};
      if(@$total_inputs != @$completed_strings){
        print "Total inputs don't match the number of completed strings can't run\n" 
            if($verbose);
        return 0;
      }else{
         my %hash;
         foreach my $input_string(@$total_inputs){
           $hash{$input_string} = 1;
         }
         my @input_list = keys(%hash);
         %hash = ();
         foreach my $input_string(@$completed_strings){
           $hash{$input_string->other->name} = 1;
         }
         my @completed_list = keys(%hash);
         my $input_set = ReseqTrack::Tools::Intersection->new(
           -list => \@input_list,
             );
         my $completed_set = ReseqTrack::Tools::Intersection->new(
           -list => \@completed_list,
             );
         my $only_input = $input_set->not($completed_set);
         if(@{$only_input->list} >= 1){
           print "There are to many ids which are only on the input set can'trun\n"
               if($verbose);
           return 0;
         }
      }
    }
  }
  print $input." can run ".$workflow->goal_event->name." ".$event->name."\n" if($verbose);
  return 1;
}


=head2 create_job_object

  Arg [1]   : ReseqTrack::Event
  Arg [2]   : string, input string
  Arg [3]   : ReseqTrack::DBSQL::DBAdaptor
  Arg [4]   : int, number of output directories
  Arg [5]   : int, max number of retries
  Function  : create ReseqTrack::Job object and store in database
  Returntype: ReseqTrack::Job
  Exceptions: n/a
  Example   : my $job = create_job_object($event, "filename", $db, 10, 4);

=cut



sub create_job_object{
  my ($event, $input_string, $db, $num_output_dirs, $retry_max) = @_;
  my $ja = $db->get_JobAdaptor;
  my $job = $ja->fetch_by_name_and_input_string($event->name, 
                                                $input_string);
  my $update = 0;
  if($job){
    if($job->retry_count && $job->retry_count >= $retry_max){
      #print STDERR $job->input_string." ".$job->event->name." has been retried too "
      #    "many times\n";
      return undef;
    }
    my $new_count = $job->retry_count;
    $new_count++;
    $job->retry_count($new_count);
    $update = 1;
  }else{
    $job = ReseqTrack::Job->new
        (
         -input_string => $input_string,
         -event => $event,
         -retry_count => 0,
         -current_status => 'CREATED',
        );
  }
  my $stem = create_job_output_stem($job->input_string, $event, $num_output_dirs);
  $job->output_file($stem.".out") unless($job->output_file);
  $job->farm_log_file($stem.".log") unless($job->farm_log_file);
  if($update){
    $ja->update($job);
    $job->current_status('CREATED');
    $ja->set_status($job);
  }else{
    $ja->store($job);
  }
  return $job;
}


=head2 create_job_output_stem

  Arg [1]   : ReseqTrack::Event
  Arg [2]   : int, number of output directories
  Function  : create output file paths
  Returntype: string, path and stem to filename
  Exceptions: throws if fails to create appropriate directory
  Example   : my $stem = create_job_output_stem($event, 10);

=cut


sub create_job_output_stem{
  my ($input_string, $event, $num_output_dirs) = @_;
  if(!$num_output_dirs){
    $num_output_dirs = 10;
  }
  my $name = basename($input_string);
  my $num = int(rand($num_output_dirs));
  my $dir = $event->output_path . "/$num/";
  if(! -e $dir){
    eval{
      mkpath($dir);
    };
    if($@){
      throw("Failed to create ".$dir." $@");
    }
  }
  my $ident = int(rand(10000));
  my $stem = $dir."/".$name."_".$event->name."_".$ident;
  $stem =~ s/\/\//\//;
  return $stem;
}


=head2 create_submission_cmds

  Arg [1]   : arrayref of ReseqTrack::Job objects
  Arg [2]   : ReseqTrack::DBSQL::DBAdaptor
  Arg [3]   : ReseqTrack::Tools::BatchSubmission object
  Arg [4]   : path to the runner script, it should point to 
    ReseqTrack/scripts/event/runner.pl
  Arg [5]   : ReseqTrack::Event object
  Function  : create batch submission cmds for the Job objects
  Returntype: hashref, keyed on command values are an arrayref of associated job 
  objects
  Exceptions: n/a
  Example   : my $cmds = create_submission_cmds($jobs, $db, $lsf, $runner, $event);

=cut


sub create_submission_cmds{
  my ($jobs, $db, $batch_object, $runner, $event) = @_;
  my $wrapper = "perl ".$runner." -dbhost ".$db->dbc->host." -dbuser ".$db->dbc->username
            ." -dbpass ".$db->dbc->password." -dbport ".$db->dbc->port
            ." -dbname ".$db->dbc->dbname." -name ".$event->name;

  if ($event->runner_options) {
      $wrapper .= " ".$event->runner_options;
  }

  my $job_count = 0;
  my @submission_jobs;

  my $ja = $db->get_JobAdaptor;
  my %cmds;
  
  my $submission_wrapper = $wrapper;
  my %jobs;
  JOB:
  foreach my $job(@$jobs){
      if(!$job->dbID){
          $job->current_status('CREATED');
          $ja->store($job);
      }
      unless($jobs{$job->dbID}){
          $jobs{$job->dbID} = 1;
      }
      else{
          throw("Duplicate version of ".$job->dbID);
      }
      $submission_wrapper .= " ".$job->dbID;

      push(@submission_jobs, $job);
      $job_count++;

      if (scalar @submission_jobs >= $event->max_array_size || $job_count >= @$jobs) {
          foreach my $submission_job (grep {$_ != $job} @submission_jobs) {
              $submission_job->farm_log_file($job->farm_log_file);
              $ja->update($submission_job);
          }
          my $submission_array_size = $event->max_array_size > 0 ? @submission_jobs : 0;
          my $cmd = $batch_object->construct_command_line($submission_wrapper, 
                                                          $event->farm_options,
                                                          $job->farm_log_file,
                                                          $event->name,
                                                          $submission_array_size,
                                                          $event->job_slot_limit);
          push(@{$cmds{$cmd}}, @submission_jobs);

          @submission_jobs = ();
          $submission_wrapper = $wrapper;
      }
  }

  return(\%cmds);
}



=head2 has_to_many_jobs

  Arg [1]   : int, maximum number of jobs
  Arg [2]   : ReseqTrack::Tools::BatchSubmission
  Arg [3]   : int, number of seconds to sleep for
  Arg [4]   : string, username to consider, defaults to $ENV{'USER'}
  Function  : checks if there are to many jobs in the submission system
  Returntype: 0
  Exceptions: n/a
  Example   : 

=cut


sub has_to_many_jobs{
  my ($max, $batch_object, $sleep, $user) = @_;
  $user = $ENV{'USER'} unless($user);
  my $total = $batch_object->get_total_job_number($user);
  if($total < $max){
    return 0;
  }else{
    print "Have ".$total." jobs compared to ".$max." going to sleep for ".$sleep.
        " seconds\n";
    sleep($sleep);
    return has_to_many_jobs($max, $batch_object, $sleep, $user);
  }
}


=head2 check_if_done

  Arg [1]   : ReseqTrack::DBSQL::DBAdaptor
  Arg [2]   : int, maximum number of times a job can be retried
  Function  : state if a pipeline is finished. This works on the assumption that 
  a finished pipeline either has no jobs left or only jobs which can no longer be
  retried
  Returntype: 0/1 
  Exceptions: n/a
  Example   : if(check_if_done($db, 4));

=cut


sub check_if_done{
  my ($db, $max_retry) = @_;
  my $ja = $db->get_JobAdaptor;
  my $jobs = $ja->fetch_all;
 JOB:foreach my $job(@$jobs){
   my $status = $job->current_status(undef, 1);
   if($status eq 'SUCCESFUL' || $status eq 'FAIL_NO_RETRY'){
      next JOB;
   }elsif($status eq 'FAILED' || $status eq 'AWOL'){ 
     if($job->retry_count < $max_retry){
       return 0;
     }else{
       next JOB;
     }
   }else{
     return 0;
   }
  }
  return 1;
}





1;
