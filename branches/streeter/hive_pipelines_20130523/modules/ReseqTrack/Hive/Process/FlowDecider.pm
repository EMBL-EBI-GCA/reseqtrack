
package ReseqTrack::Hive::Process::FlowDecider;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use List::Util qw(sum);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

#    my $flows = [];
#    if ($self->param_is_defined('count_files') && $self->param('count_files')) {
#      my $num_files = $self->count_param_values('files');
#
#      my $flows_if_no_files = $self->get_param_values('flows_if_no_files');
#      my $flows_if_one_file = $self->get_param_values('flows_if_one_file');
#      my $flows_if_multiple_files = $self->get_param_values('flows_if_multiple_files');
#
#      $flows = scalar $num_files >1 ? $flows_if_multiple_files
#                      : scalar $num_files ==1 ? $flows_if_one_file
#                      : $flows_if_no_files;
#
#      if ($self->param_is_defined('require_true')) {
#        my $requires = $self->param('require_true');
#        $flows = [grep {$requires->{$_} // 1} @$flows];
#      }
#    }
#    elsif ($self->param_is_defined('require_true')) {
#      my $requires = $self->param('require_true');
#      $flows = [grep {$requires->{$_}} keys %$requires];
#    }


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

