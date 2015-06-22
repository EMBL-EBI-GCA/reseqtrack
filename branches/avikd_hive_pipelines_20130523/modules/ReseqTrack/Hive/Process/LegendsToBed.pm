package ReseqTrack::Hive::Process::LegendsToBed;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_executable check_file_exists);
use ReseqTrack::Hive::Utils::ImputeSliceUtils qw( reference_config_to_bed_slices filter_haps_missing );
use ReseqTrack::Hive::Utils::SequenceSliceUtils qw( fai_to_slices );

=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    my $reference_config = $self->param_required('reference_config');
    my $fai = $self->param_required('fai');
    check_file_exists($reference_config);
    
    $self->param_required('shapeit_haps');
    my $haps = $self->file_param_to_flat_array('shapeit_haps');
    
    throw("expecting single haps") if(scalar @$haps >1);
    
    my $hap = @$haps[0];
    
    my $remove_haps_missing = $self->param_is_defined('remove_haps_missing') ? $self->param('remove_haps_missing') : undef ;
    
    my $max_base = $self->param_is_defined('max_base') ? $self->param('max_base') : 5000000 ;
    
    my $SQ_start = $self->param_required('SQ_start');
    my $SQ_end = $self->param_required('SQ_end');
    my $bp_start = $self->param_required('bp_start');
    my $bp_end = $self->param_required('bp_end');
    
    throw("expecting single chromosome") unless($SQ_start eq $SQ_end);
    
    my $chrX_string = $self->param_is_defined('chrX_string') ? $self->param('chrX_string') : undef;
    
    my $chrX_PAR1_start = $self->param_is_defined('chrX_PAR1_start') ? $self->param('chrX_PAR1_start') : 1 ;
    my $chrX_PAR1_end = $self->param_is_defined('chrX_PAR1_end') ? $self->param('chrX_PAR1_end') : 2699500 ;
    my $chrX_PAR2_start = $self->param_is_defined('chrX_PAR2_start') ? $self->param('chrX_PAR2_start') : 154933000 ;
    my $chrX_PAR2_end = $self->param_is_defined('chrX_PAR2_end') ? $self->param('chrX_PAR2_end') : 155270560 ;
    
    my $gunzip = $self->param_is_defined('gunzip') ? $self->param('gunzip') : "gunzip" ;
    
    my $output_dir = $self->output_dir;
    my $job_name = $self->job_name;
    
    my $slices = fai_to_slices(
          fai => $fai,
          SQ_start => $SQ_start, SQ_end => $SQ_end,
          bp_start => $bp_start, bp_end => $bp_end,
          );
    
    throw("expecting single slice, got :", scalar @$slices) if(scalar @$slices > 1);
    
    my $slice = $$slices[0];
    
    check_directory_exists($output_dir);
    
    my $output_file = $output_dir. '/'. $SQ_start;
    $output_file .= '.'. $bp_start .'-'. $bp_end if(!$slice->is_whole_SQ); 
    $output_file .='.bed';
    
    open my $OUT, "> $output_file";   
  
    my $bed_slices = reference_config_to_bed_slices(
          reference_config=> $reference_config,
          SQ_start => $SQ_start,
          bp_start => $bp_start,
          bp_end => $bp_end,
          max_base => $max_base,
          gunzip => $gunzip,
          chrX_PAR1_start => $chrX_PAR1_start, 
          chrX_PAR1_end => $chrX_PAR1_end, 
          chrX_PAR2_start => $chrX_PAR2_start, 
          chrX_PAR2_end => $chrX_PAR2_end,
          chrX_string => $chrX_string,
    );
    
    foreach my $bed_slice (@$bed_slices) {
      
      if( $remove_haps_missing == 1) {
      
        $bed_slice = filter_haps_missing(
                  haps => $hap,
                  slice => $bed_slice,
                  SQ_start => $SQ_start,
                );
      }
      
      print $OUT join "\n", @$bed_slice ;
    }
    
    close $OUT;
    
    $self->output_param('legends_bed' , $output_file);  
}

1;