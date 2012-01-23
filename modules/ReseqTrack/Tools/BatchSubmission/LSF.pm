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
use List::Util qw(first);

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
  Arg [4]   : string, file path to use for LSF output
  Arg [5]   : string to use as job name if not already part of options
  Arg [6]   : int, size of lsf job array (0 implies it is not an array)
  Function  : construct bsub command
  Returntype: string, bsub command
  Exceptions: n/a
  Example   : my $bsub_cmd = $lsf->construct_command_line($to_run, $options, $farm_log_file, $name, $array_size);
                                                            

=cut



sub construct_command_line{
  my ($self, $to_run, $options, $farm_log_file, $name, $array_size, $job_slot_limit) = @_;

  $to_run = $self->cmd unless($to_run);
  $options = $self->options unless($options);

  $options =~ s/\-J\s+\S+//;
  $options .= ' -J '.$name;
  if ($array_size >0) {
      $options .= '[1-'.$array_size.']';
      if ($job_slot_limit) {
          $options .= '%'.$job_slot_limit;
      }
  }

  my $cmd = $self->program." ".$options." -o ".$farm_log_file." ".$to_run;
  #print "LSF cmd ".$cmd."\n";
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


sub refresh_job_info{
    my ($self, $submission_id, $submission_index) = @_;

    my $cmd = 'bjobs -l ' . $submission_id;
    if ($submission_index) {
        $cmd .= "[" . $submission_index . "]";
    }

  open(CMD, $cmd." |") or throw("ReseqTrack::Tools::BatchSubmission::LSF ".
                                "Failed to run $cmd $!");

  my @job_info;

  $cmd .= "\n";
  push(@job_info, $cmd);

  while (my $line = <CMD>) {
      push(@job_info, $line);
  }
  close CMD;

  $self->job_info($submission_id, $submission_index, undef, \@job_info);

  return \@job_info;
}

=head2 job_info

  Arg [1]   : ReseqTrack::Tools::BatchSubmission::LSF
  Arg [2]   : int, submission_id
  Arg [3]   : int, submission_index
  Arg [4]   : boolean, refresh, to refetch job_info using method refresh_job_info
  Arg [5]   : arrayref, optional, output from 'bjobs -l id[index]'
  Function  : accessor method for output from 'bjobs -l id[index]'
  Returntype: arrayref of strings containing all of the output from 'bjobs -l id[index]'
  Exceptions: throws if fails to run command
  Example   : ReseqTrack::Tools::BatchSubmission::LSF->job_info($id, $index)

=cut

sub job_info{
  my ($self, $submission_id, $submission_index, $refresh, $job_info) = @_;

  if ($job_info) {
    $self->{'job_info'}->{$submission_id}->{$submission_index} = $job_info;
  }
  elsif (! $self->{'job_info'}->{$submission_id}->{$submission_index} || $refresh) {
    $self->refresh_job_info($submission_id, $submission_index);
  }
  return $self->{'job_info'}->{$submission_id}->{$submission_index};
}

=head2 memory_usage

  Arg [1]   : ReseqTrack::Tools::BatchSubmission::LSF
  Arg [2]   : int, submission_id
  Arg [3]   : int, submission_index
  Arg [4]   : boolean, refresh, to refetch job_info using method refresh_job_info
  Function  : get memory_usage and swap_usage from 'bjobs -l id[index]'
  Returntype: array of two integers, memory_usage and swap_usage
  Exceptions: throws if fails to run command
  Example   : ($mem, $swap) = ReseqTrack::Tools::BatchSubmission::LSF->memory_usage($id, $index)

=cut

sub memory_usage{
  my ($self, $submission_id, $submission_index, $refresh) = @_;

  my $job_info = $self->job_info($submission_id, $submission_index, $refresh);
  my $mem_line = first { /^ *MEM:.*; *SWAP:.*;/ } @$job_info;
  throw("could not find memory usage") if (!$mem_line);

  $mem_line =~ /^ *MEM: *(.*); *SWAP: *(.*);/;
  my ($memory, $swap) = ($1, $2);
  throw("could not find memory usage") if (!$memory || !$swap);

  foreach ($memory, $swap) {
    /(\d+) *(\w+)/;
    my $value = $1;
    my $units = $2;
    $value = $units eq 'Gbytes' ? $value * 1024
          : $units eq 'Mbytes' ? $value
          : $units eq 'Kbytes' ? int ($value/1024 + 0.5)
          : 0;
    $_ = $value;
  }

  return ($memory, $swap);
}


1;
