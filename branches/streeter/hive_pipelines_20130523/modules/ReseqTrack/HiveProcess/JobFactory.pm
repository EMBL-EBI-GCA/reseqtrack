
package ReseqTrack::HiveProcess::JobFactory;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $values = $self->get_param_values('data_value');

    foreach my $value (@$values) {
      $self->prepare_child_output_id($value, {'data_value' => $value});
    }

}

1;

