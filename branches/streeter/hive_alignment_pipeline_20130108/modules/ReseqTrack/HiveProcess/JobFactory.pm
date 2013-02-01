
package ReseqTrack::HiveProcess::JobFactory;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $values = $self->param('values') || die "'value' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";
    my $branch_label = $self->param('branch_label');

    $values = ref($values) eq 'ARRAY' ? $values : [$values];
    foreach my $value (@$values) {
      my $new_branch_label = $branch_label ? "$branch_label.$value" : $value;
      $self->output_child_branches('value' => $value, 'label' => $new_branch_label, 'output_dir' => "$output_dir/$value");
    }

}

1;

