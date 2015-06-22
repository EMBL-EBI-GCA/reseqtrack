package ReseqTrack::Hive::Process::RunMergePhasedVcf;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw( check_directory_exists check_file_exists );
use File::Copy qw(move);
use IPC::System::Simple qw(system);
use ReseqTrack::Tools::RunMergePhasedVcf;


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    $self->param_required('vcf');
    $self->param_required('SQ_start');
    $self->param_required('SQ_end');
    
    my $bgzip = $self->param_is_defined('bgzip') ? $self->param('bgzip') : 'bgzip';
    my $tabix = $self->param_is_defined('tabix') ? $self->param('tabix') : 'tabix';
    my $run_tabix = $self->param_is_defined('run_tabix') && $self->param('run_tabix') ? 1 : 0;
    
    my $vcfs = $self->file_param_to_flat_array('vcf');
    my $SQ_start = $self->param_to_flat_array('SQ_start');
    my $SQ_end = $self->param_to_flat_array('SQ_end');
    
    foreach my $vcf_path (grep {defined $_} @$vcfs) {
      check_file_exists($vcf_path);
    }
    
    my $output_dir = $self->output_dir;
    my $job_name = $self->job_name;
    my $output_file = "$output_dir/$job_name.vcf.gz";
    check_directory_exists($output_dir);
    
    if (scalar @$vcfs == 1) {
      move($$vcfs[0], $output_file) or throw("could not move $$vcfs[0] to $output_file $!");
      if ($run_tabix) {
        system("$tabix -p vcf $output_file") ==0 or throw("tabix failed $!");
      }
    }
    elsif (scalar @$vcfs > 1) {
      
      my $vcf_list;
      my $chrom = $SQ_start->[0];
      
      foreach my $i (0..$#{$vcfs}) {
        throw("expecting single chromosome, got", $SQ_start->[$i]," and ",$SQ_end->[$i]) unless($SQ_start->[$i] eq $SQ_end->[$i]);
        throw("got a different chromosome,",$SQ_start->[$i]," expecting $chrom") unless($SQ_start->[$i] eq $chrom);
        
        push @{$vcf_list}, $vcfs->[$i];
      }
      
      my $merge_object = ReseqTrack::Tools::RunMergePhasedVcf->new(
        -input_files  => $vcf_list,
        -program      => $self->param_is_defined('program_file') ? $self->param('program_file') : undef,
        -working_dir  => $self->output_dir,
        -job_name     => $self->job_name,
        -chrom        => $chrom,
        -bgzip        => $self->param_is_defined('bgzip') ? $self->param('bgzip') : undef,
        -tabix        => $self->param_is_defined('tabix') ? $self->param('tabix') : undef,
        -create_index => $run_tabix,
      );

      $self->run_program($merge_object);
      my $output_files = $merge_object->output_vcf_files;
      $output_file = $output_files->[0];
      
#       my $output_vcf_files = $merge_object->output_vcf_files;
#       move($$output_vcf_files[0], $output_file) or throw("could not move $$output_vcf_files[0] to $output_file $!");
        
    }
    else {
      throw("no vcf file to merge");
    }
    
    
    
    $self->output_param('vcf' => $output_file);
    if ($run_tabix) {
      $self->output_param('tbi' => "$output_file.tbi");
    }
    
}

1;