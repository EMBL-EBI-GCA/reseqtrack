
package ReseqTrack::HiveProcess::MergeVcf;

use strict;

use base ('ReseqTrack::HiveProcess::BranchableProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_executable check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $vcfs = $self->get_param_array('vcf');
    my $bgzip = $self->param('bgzip') || 'bgzip';

    throw("no vcf files") if !scalar @$vcfs;
    foreach my $vcf_path (@$vcfs) {
      check_file_exists($vcf_path);
    }

    my $output_dir = $self->output_dir;
    my $job_name = $self->job_name;
    my $output_file = "$output_dir/$job_name.vcf.gz";
    check_directory_exists($output_dir);
    check_executable($bgzip);
    open my $OUT, "| $bgzip -c > $output_file";
    my $first_vcf = 1;
    $self->data_dbc->disconnect_when_inactive(1);
    foreach my $vcf_path (@$vcfs) {

      ##############
      #This is a temporary hack:
      #Needs to be fixed so I get this from database.
      use File::Basename qw(basename);
      my ($region_start, $region_end) = basename($vcf_path) =~ /\.(\d+)-(\d+)\./;
      ##############

      my $IN;
      if ($vcf_path =~ /\.b?gz(?:ip)?$/){
        open $IN, "$bgzip -cd $vcf_path |" or throw("cannot open $vcf_path: $!");
      }
      else {
        open $IN, '<', $vcf_path or throw("cannot open $vcf_path: $!");
      }
      LINE:
      while (my $line = <$IN>) {
        if ($line =~ /^#/) {
          print $OUT $line if $first_vcf;
          next LINE;
        }
        my ($pos) = $line =~ /\t(\d+)/;
        next LINE if $pos < $region_start;
        next LINE if $pos > $region_end;
        print $OUT $line;
      }
      close $IN;
      $first_vcf = 0;
    }
    close $OUT;
    $self->data_dbc->disconnect_when_inactive(0);

    $self->output_this_branch('vcf' => $output_file);

}


1;

