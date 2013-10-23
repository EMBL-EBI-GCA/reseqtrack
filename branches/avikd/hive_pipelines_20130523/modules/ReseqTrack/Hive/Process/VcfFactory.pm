
package ReseqTrack::Hive::Process::VcfFactory;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use POSIX qw(ceil);
use ReseqTrack::Tools::Exception qw(throw);
use File::Basename qw(fileparse);

=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    
    $self->param_required('vcf_list');
    
      my $vcfs = $self->file_param_to_flat_array('vcf_list');
      
      my $count=0;
      
      foreach my $vcf_path (@{$vcfs}) {
               
       
        my $label = fileparse($vcf_path, qw( .vcf .vcf.gz )); 
       
        $self->prepare_factory_output_id("$label", {
                                                    'vcf' => $vcf_path,
                                                    'fan_index' => $count,});
        $count++;
     
        }
      
    

}

1;

