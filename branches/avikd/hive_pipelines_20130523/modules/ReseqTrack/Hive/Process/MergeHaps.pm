
package ReseqTrack::Hive::Process::MergeHaps;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_executable check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    $self->param_required('hap');
    $self->param_required('bp_start');
    $self->param_required('bp_end');
    my $zcat = $self->param_is_defined('zcat') ? $self->param('zcat') : 'zcat';
  

    my $haps = $self->file_param_to_flat_array('hap');
    my $bp_start = $self->param_to_flat_array('bp_start');
    my $bp_end = $self->param_to_flat_array('bp_end');
    foreach my $hap_path (grep {defined $_} @$haps) {
      check_file_exists($hap_path);
    }

    my $output_dir = $self->output_dir;
    my $job_name = $self->job_name;
    my $output_file = "$output_dir/$job_name.haps";
    check_directory_exists($output_dir);
    check_executable($zcat);
    
    open my $OUT, " > $output_file";
    my $first_hap = 1;
    $self->dbc->disconnect_when_inactive(1);
    HAPS:
    foreach my $i (0..$#{$haps}) {
      my $hap_path = $haps->[$i];
      
      next HAPS if ! defined $hap_path;
      throw("missing bp_start for $hap_path") if ! defined $bp_start->[$i];
      throw("missing bp_end for $hap_path") if ! defined $bp_end->[$i];

      my $IN;

      open $IN, '<', $hap_path or throw("cannot open $hap_path: $!");
      
      LINE:
      while (my $line = <$IN>) {
        if ($line =~ /^#/) {
          print $OUT $line if $first_hap;
          next LINE;
        }
        my ($pos) = $line =~ /\S+\s+(\d+)/; 
        next LINE if $pos < $bp_start->[$i];
        next LINE if $pos > $bp_end->[$i];
        print $OUT $line;
      }
      close $IN;
      $first_hap = 0;
    }
    close $OUT;

    $self->dbc->disconnect_when_inactive(0);

    $self->output_param('hap' => $output_file);

}


1;

