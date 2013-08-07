
package ReseqTrack::Hive::Process::MergeVcf;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_executable check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    $self->param_required('vcf');
    $self->param_required('bp_start');
    $self->param_required('bp_end');
    my $bgzip = $self->param_is_defined('bgzip') ? $self->param('bgzip') : 'bgzip';
    my $tabix = $self->param_is_defined('tabix') ? $self->param('tabix') : 'tabix';
    my $run_tabix = $self->param_is_defined('run_tabix') && $self->param('run_tabix') ? 1 : 0;

    my $vcfs = $self->file_param_to_flat_array('vcf');
    my $bp_start = $self->param_to_flat_array('bp_start');
    my $bp_end = $self->param_to_flat_array('bp_end');
    foreach my $vcf_path (grep {defined $_} @$vcfs) {
      check_file_exists($vcf_path);
    }

    my $output_dir = $self->output_dir;
    my $job_name = $self->job_name;
    my $output_file = "$output_dir/$job_name.vcf.gz";
    check_directory_exists($output_dir);
    check_executable($bgzip);
    if ($run_tabix) {
      check_executable($tabix);
    }
    open my $OUT, "| $bgzip -c > $output_file";
    my $first_vcf = 1;
    $self->dbc->disconnect_when_inactive(1);
    VCF:
    foreach my $i (0..$#{$vcfs}) {
      my $vcf_path = $vcfs->[$i];
      next VCF if ! defined $vcf_path;
      throw("missing bp_start for $vcf_path") if ! defined $bp_start->[$i];
      throw("missing bp_end for $vcf_path") if ! defined $bp_end->[$i];

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
        next LINE if $pos < $bp_start->[$i];
        next LINE if $pos > $bp_end->[$i];
        print $OUT $line;
      }
      close $IN;
      $first_vcf = 0;
    }
    close $OUT;

    if ($run_tabix) {
      system("$tabix -p vcf $output_file") ==0 or throw("tabix failed $!");
    }

    $self->dbc->disconnect_when_inactive(0);

    $self->output_param('vcf' => $output_file);
    if ($run_tabix) {
      $self->output_param('tbi' => "$output_file.tbi");
    }

}


1;

