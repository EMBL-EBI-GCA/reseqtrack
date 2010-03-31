=pod

=head1 NAME

ReseqTrack::Tools::BatchSubmission::LSF

=head1 SYNOPSIS

This object provides functionality to create LSF bsub command lines and also
check on job status

=head1 Example

my $lsf = ReseqTrack::Tools::BatchSubmission::LSF->new(
                                                       -program => "bsub",
                                                       -max_job_number => 10000,
                                                       -sleep => 360,
                                                      );
my $bsub_cmd = $lsf->construct_comment_line($cmd_to_run, "-q production -J job_name",
                                            $job_object)l

=cut


package ReseqTrack::Tools::BatchSubmission::LSF;

use strict;
use warnings;
use vars qw(@ISA);


use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::BatchSubmission;

@ISA = qw(ReseqTrack::Tools::BatchSubmission);



=head2 new

  Arg [1]   : ReseqTrack::Tools::BatchSubmission::LSF
  Function  : create LSF object. By default it sets the program to be bsub and the
  user to be $ENV{'USER'}
  Returntype: ReseqTrack::Tools::BatchSubmission::LSF
  Exceptions: 
  Example   : 

=cut


sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);

  ####default setting####
  $self->program("bsub") unless($self->program);
  $self->user($ENV{'USER'}) unless($self->user);
  ######################

  return $self;
}


=head2 construct_command_line

  Arg [1]   : ReseqTrack::Tools::BatchSubmission::LSF
  Arg [2]   : string, command line that should be submitted 
  Arg [3]   : string, commandline options for bsub
  Arg [4]   : ReseqTrack::Job object associated with job
  Arg [5]   : string to use as job name if not already part of options
  Function  : construct bsub command
  Returntype: string, bsub command
  Exceptions: n/a
  Example   : my $bsub_cmd = $lsf->construct_command_line($to_run, $options, $job);
                                                            

=cut



sub construct_command_line{
  my ($self, $to_run, $options, $job, $name) = @_;
  
  $to_run = $self->cmd unless($to_run);
  $options = $self->options unless($options);
  if($name){
    $options .= " -J ".$name unless($options =~ /\-J\s+\S+/);
  }
  my $cmd = $self->program." ".$options." -o ".$job->stdout_file." -e ".$job->stderr_file.
      " ".$to_run;
  return $cmd;
}


=head2 get_total_job_number

  Arg [1]   : ReseqTrack::Tools::BatchSubmission::LSF
  Arg [2]   : string, username
  Function  : return total number of jobs using busers command
  Returntype: int
  Exceptions: throws if failed to run busers 
  Example   : 

=cut



sub get_total_job_number{
  my ($self, $user) = @_;
  $user = $self->user unless($user);
  my $cmd = "busers $user";
  open(FH, $cmd." | ")  or throw("ReseqTrack::Tools::BatchSubmission::LSF ".
                                 "Failed to open ".$cmd." $!");
  my $total;
  while(<FH>){
    next if(/USER\/GROUP/);
    chomp;
    my @values = split;
    $total = $values[4];
  }
  close(FH);
  return $total;
}



=head2 check_for_awol_jobs

  Arg [1]   : ReseqTrack::Tools::BatchSubmission::LSF
  Arg [2]   : arrayref of ReseqTrack::Job objects
  Arg [3]   : username
  Function  : ensures all current pipeline jobs still exist in LSF
  Returntype: n/a
  Exceptions: 
  Example   : 

=cut



sub check_for_awol_jobs{
  my ($self, $jobs, $user) = @_;
  $user = $self->user unless($user);
  my $ja = $self->db->get_JobAdaptor;
  $jobs = $ja->fetch_all unless($jobs && @$jobs >= 1);
  my $submission_hash = run_bjobs($user);
  foreach my $job(@$jobs){
    next if($job->current_status =~ /^FAIL/);
    next if($job->current_status =~ /AWOL/);
    next if($job->current_status eq 'CREATED');
    next if(!$job->submission_id || $job->submission_id == 0);
    unless($submission_hash->{$job->submission_id}){
      $job->current_status("AWOL");
      $ja->set_status($job);
    }
  }
} 


=head2 run_bjobs

  Arg [1]   : ReseqTrack::Tools::BatchSubmission::LSF
  Arg [2]   : string, username
  Function  : run bjobs cmd
  Returntype: hash keyed on submission id
  Exceptions: throws if fails to run command
  Example   : 

=cut



sub run_bjobs{
  my ($self, $user) = @_;
  $user = $self->user unless($user);
  my $cmd = "bjobs -w -u $user";
  open(CMD, $cmd." |") or throw("ReseqTrack::Tools::BatchSubmission::LSF ".
                                "Failed to run ".$cmd);
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
  close(CMD);
  return \%hash;
}

1;
