
package ReseqTrack::HiveProcess::JobFactory;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $values = $self->get_param_array('data_values');

    foreach my $value (@$values) {
      $self->output_child_branches('data_value' => $value, 'label' => $value);
    }

}

1;

