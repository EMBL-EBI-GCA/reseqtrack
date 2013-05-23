
package ReseqTrack::HiveProcess::FlowDecider;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $files = $self->get_param_array('file');
    my $flows_if_no_files = $self->param('flows_if_no_files');
    my $flows_if_one_file = $self->param('flows_if_one_file');
    my $flows_if_multiple_files = $self->param('flows_if_multiple_files');

    my $flows = scalar @$files >1 ? $flows_if_multiple_files
                    : scalar @$files >0 ? $flows_if_one_file
                    : $flows_if_no_files;

    $self->flows_this_branch($flows);

}


1;

