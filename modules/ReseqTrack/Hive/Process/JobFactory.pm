
package ReseqTrack::Hive::Process::JobFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use File::Basename qw(fileparse);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

   # my $values = $self->param_required('factory_value');
    my $values = $self->file_param_to_flat_array('factory_value');

    foreach my $value (@$values) {
      if($value=~ /^\//) {
        my $label = fileparse($value); 
        $self->prepare_factory_output_id($label, {'factory_value' => $value});
        
      }
      else {       
        $self->prepare_factory_output_id($value, {'factory_value' => $value});
      }
    }

}

1;

