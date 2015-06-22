
package ReseqTrack::Hive::Process::JobFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $values = $self->param_required('factory_value');

    foreach my $value (@$values) {
      $self->prepare_factory_output_id({'factory_value' => $value});
    }

}

1;

