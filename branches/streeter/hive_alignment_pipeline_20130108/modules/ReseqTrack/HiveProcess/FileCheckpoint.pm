
package ReseqTrack::HiveProcess::FileCheckpoint;

use strict;
use warnings;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use File::Basename qw( fileparse );


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $files = $self->param('file') || die "'file' is an obligatory parameter";
    my $move_to_branch_directory = $self->param('move_to_branch_directory');
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";
    my $branch_label = $self->param('branch_label') || die "'branch_label' is an obligatory parameter";
    my $set_never_delete = $self->param('set_never_delete');
    my $reset_name = $self->param('reset_name');
    my $reset_dir = $self->param('reset_dir');
    my $suffix = $self->param('suffix');

    $files = ref($files) eq 'ARRAY' $files : [$files];
    foreach my $file (@$files) {
      my ($name, $dir) = fileparse($file);
      my ($new_name, $new_dir) = ($name, $dir);
      if ($reset_name) {
        $new_name = "$branch_label.$suffix";
      }
      if ($reset_dir) {
        $new_dir = $output_dir;
      }
      my $new_path = $new_dir.$new_name;
    }

}


1;

