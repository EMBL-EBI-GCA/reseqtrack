
package ReseqTrack::HiveProcess::SplitFastq;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_file_exists);
use ReseqTrack::Tools::Exception qw(throw);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;

    my $run_id = $self->param('run_id') || die "'run_id' is an obligatory parameter";
    my $type_fastq = $self->param('type_fastq') || die "'type_fastq' is an obligatory parameter";
    my $type_split_fastq = $self->param('type_split_fastq') || die "'type_split_fastq' is an obligatory parameter";
    my $max_bases = $self->param('max_bases') || die "'max_bases' is an obligatory parameter";
    my $output_dir = $self->param('output_dir') || die "'output_dir' is an obligatory parameter";
    my $program_file = $self->param('program_file');
    my $directory_layout = $self->param('directory_layout');

    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$self->param('reseqtrack_db')});
    my $ca = $db->get_CollectionAdaptor;
    my $collection = $ca->fetch_by_name_and_type($run_id, $type_fastq);
    throw("Failed to find a collection for $run_id $type_fastq") if(!$collection);
    my @input_fastq = map {$_->name} @{$collection->others};

    my $rmi_a = $db->get_RunMetaInfoAdaptor;
    my $run_meta_info = $rmi_a->fetch_by_run_id($run_id);
    throw("Failed to find run_meta_info for $run_id") if (!$run_meta_info);

    foreach my $fastq (@input_fastq) {
      check_file_exists($fastq);
    }

    check_directory_exists($output_dir);

    foreach (1..5) {
      my $mate1 = "$output_dir/${run_id}_1.$_.fastq";
      my $mate2 = "$output_dir/${run_id}_2.$_.fastq";
      system("touch $mate1");
      system("touch $mate2");
      $self->output_child_branches($type_split_fastq => [$mate1, $mate2], 'label' => "$run_id.$_", 'output_dir' => "$output_dir/$_");
    }

}

1;

