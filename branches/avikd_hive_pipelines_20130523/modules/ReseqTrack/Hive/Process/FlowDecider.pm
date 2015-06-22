=head1 NAME

ReseqTrack::Hive::Process::FlowDecider;

=head1 SYNOPSIS

    Hive process to control flow of the pipeline.

    Best illustrated with an example.  The configuration parameters for decide_merge_chunks from the Alignment pipeline:

          -parameters => {
              realign_knowns_only => $self->o('realign_knowns_only'),
              recalibrate_level => $self->o('recalibrate_level'),
              files => '#bam#',
              require_file_count => {
                        1 => '1+',
                        2 => '2+',
                        3 => '1+',
                      },
              require_true => {
                  1 => '#expr($realign_knowns_only || $recalibrate_level==1)expr#',
                  2 => '#expr($realign_knowns_only || $recalibrate_level==1)expr#',
                  3 => '#expr(!$realign_knowns_only && $recalibrate_level!=1)expr#',
              }
          },

    This module looks at the parameter named 'files'.  Here, files is set equal to 'bam' i.e. it counts the number of bam files passed to this job.
    The parameters 'realign_knowns_only' and 'recalibrate_level' are used as boolean values to control the flow

    This job may flow down flows 1, 2, and/or 3.

    It will flow down #1 if there is 1 or more bam file and if ($realign_knowns_only || $recalibrate_level==1)
    It will flow down #2 if there is 2 or more bam file and if ($realign_knowns_only || $recalibrate_level==1)
    It will flow down #3 if there is 1 or more bam file and if (!$realign_knowns_only || $recalibrate_level==1)

    When using this module, 'require_true' and 'require_file_count' are both optional (but you will want to define at least one of them)


=cut

package ReseqTrack::Hive::Process::FlowDecider;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use List::Util qw(sum);

sub run {
    my $self = shift @_;


    my $require_true_hash = $self->param_is_defined('require_true') ? $self->param('require_true') : {};
    my $require_count_hash = {};
    my $num_files;
    if ($self->param_is_defined('require_file_count')) {
      $require_count_hash = $self->param('require_file_count');
      $num_files = $self->count_param_values('files');
    }

    my %possible_flows;
    foreach my $flow (keys %$require_true_hash, keys %$require_count_hash) {
      $possible_flows{$flow} = 1;
    }

    my @output_flows;
    FLOW:
    foreach my $flow (keys %possible_flows) {
      if (defined $require_true_hash->{$flow}) {
        next FLOW if !$require_true_hash->{$flow};
      }
      if (defined $require_count_hash->{$flow}) {
        my $condition = $require_count_hash->{$flow};
        my ($require_num, $modifier) = $require_count_hash->{$flow} =~ /(\d+)([+-]?)/;
        throw("did not recognise condition for $flow") if !defined $require_num;
        next FLOW if !$modifier && $require_num != $num_files;
        next FLOW if $modifier eq '+' && $require_num > $num_files;
        next FLOW if $modifier eq '-' && $require_num < $num_files;
      }
      push(@output_flows, $flow);
    }
    

    $self->flows_non_factory(\@output_flows);

}


1;

