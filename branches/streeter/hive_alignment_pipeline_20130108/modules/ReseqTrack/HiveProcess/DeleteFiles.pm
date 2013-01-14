
package ReseqTrack::HiveProcess::DeleteFiles;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::FileSystemUtils qw(delete_file);
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::DBSQL::DBAdaptor;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $files = $self->param('files') || die "'files' is an obligatory parameter";
    my $root_output_dir = $self->param('root_output_dir') || die "'root_output_dir' is an obligatory parameter";

    $files = ref($files) eq 'ARRAY' ? $files : [$files];
    foreach my $file (@$files) {
      throw("will not delete file outside of the root_output_dir $file $root_output_dir") if $file !~ /^$root_output_dir/;
      delete_file($file);
    }
}

1;

