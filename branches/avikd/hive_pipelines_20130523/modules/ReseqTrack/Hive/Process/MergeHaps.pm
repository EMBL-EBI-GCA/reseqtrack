
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
    $self->param_required('samples');
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
    my $output_haps = "$output_dir/$job_name.haps";
    check_directory_exists($output_dir);
    check_executable($zcat);
    
    ### HAPS
    open my $OUT, " > $output_haps" || throw("can't write $output_haps: $!") ;
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
        my ($pos) = $line =~ /^\S+\s+\S+\s+(\d+)/; 
        next LINE if $pos < $bp_start->[$i];
        next LINE if $pos > $bp_end->[$i];
        print $OUT $line;
      }
      close $IN;
      $first_hap = 0;
    }
    close $OUT;  
      
    ### SAMPLES  
    my $samples = $self->file_param_to_flat_array('samples');
    foreach my $sample_path (grep {defined $_} @$samples) {
      check_file_exists($sample_path);
    }
    
    my $output_sample = "$output_dir/$job_name.samples";
    open my $OUT_SAMPLE, " > $output_sample" || throw("can't write $output_sample: $!");
    my $first_sample = 1;
    SAMPLE:
    foreach my $i (0..$#{$samples}) {
      my $sample_path = $samples->[$i];
      next SAMPLE if ! defined $sample_path;
      
      throw("missing bp_start for $sample_path") if ! defined $bp_start->[$i];
      throw("missing bp_end for $sample_path") if ! defined $bp_end->[$i];
      
      my $IN_sample;
      open $IN_sample, '<', $sample_path or throw("cannot open $sample_path: $!");
      
      SAMPLE_LINE:
      while (my $sample_line = <$IN_sample>) {
        print $OUT_SAMPLE $sample_line if $first_sample;
      }
      close $IN_sample;
      $first_sample = 0;
    }
    close $OUT_SAMPLE;
    
    #### IMPUTE
    my $output_impute = "$output_dir/$job_name.impute";
    
    if($self->param_is_defined('impute')) {
      my $imputes = $self->file_param_to_flat_array('impute');
      
      foreach my $impute_path (grep {defined $_} @$imputes) {
        check_file_exists($impute_path);
      }
     
      open my $OUT_IMPUTE, " > $output_impute" || throw("can't write $output_impute: $!");
      my $first_impute = 1;
      
      IMPUTE:
      foreach my $i (0..$#{$imputes}) {
        my $impute_path = $imputes->[$i];
        next IMPUTE if ! defined $impute_path;
      
        throw("missing bp_start for $impute_path") if ! defined $bp_start->[$i];
        throw("missing bp_end for $impute_path") if ! defined $bp_end->[$i];
      
        my $IN_impute;
        open $IN_impute, '<', $impute_path or throw("cannot open $impute_path: $!");
      
        IMPUTE_LINE:
        while (my $impute_line = <$IN_impute>) {
          if ($impute_line =~ /^#/) {
            print $OUT_IMPUTE $impute_line if $first_impute;
            next IMPUTE_LINE;
          }
          my ($pos) = $impute_line =~ /^\S+\s+\S+\s+(\d+)/; 
          next IMPUTE_LINE if $pos < $bp_start->[$i];
          next IMPUTE_LINE if $pos > $bp_end->[$i];
          print $OUT_IMPUTE $impute_line;
        }
        close $IN_impute;
        $first_impute = 0; 
      }
      close $OUT_IMPUTE;
   }
   
   ### SUMMARY
   my $output_summary = "$output_dir/$job_name.summary";
   
   if($self->param_is_defined('impute_summary')) {
      my $impute_summarys = $self->file_param_to_flat_array('impute_summary');
      
      foreach my $impute_summary_path (grep {defined $_} @$impute_summarys) {
        check_file_exists($impute_summary_path);
      }
      
      open my $OUT_SUMMARY, " > $output_summary" || throw("can't write $output_summary: $!");
      
      SUMMARY:
      foreach my $i (0..$#{$impute_summarys}) {
        my $impute_summary_path = $impute_summarys->[$i];
        next SUMMARY if ! defined $impute_summary_path;
        
        throw("missing bp_start for $impute_summary_path") if ! defined $bp_start->[$i];
        throw("missing bp_end for $impute_summary_path") if ! defined $bp_end->[$i];
        
        my $IN_summary;
        open $IN_summary, '<', $impute_summary_path or throw("cannot open $impute_summary_path: $!");
        
        print $OUT_SUMMARY "Summary for region:",$bp_start->[$i],"-",$bp_end->[$i],"\n\n";
        
        SUMMARY_LINE:
        while (my $summary_line = <$IN_summary>) {
           print $OUT_SUMMARY $summary_line;
        }
        close $IN_summary;
        print $OUT_SUMMARY "\n\n";
        
      }
      close $OUT_SUMMARY;
   }
   
   ### INFO
   my $output_info = "$output_dir/$job_name.info";
    
    if($self->param_is_defined('impute_info')) {
      my $impute_infos = $self->file_param_to_flat_array('impute_info');
      
      foreach my $info_path (grep {defined $_} @$impute_infos) {
        check_file_exists($info_path);
      }
     
      open my $OUT_INFO, " > $output_info" || throw("can't write $output_info: $!");
      my $first_impute = 1;
      
      INFO:
      foreach my $i (0..$#{$impute_infos}) {
        my $info_path = $impute_infos->[$i];
        next INFO if ! defined $info_path;
      
        throw("missing bp_start for $info_path") if ! defined $bp_start->[$i];
        throw("missing bp_end for $info_path") if ! defined $bp_end->[$i];
      
        my $IN_info;
        open $IN_info, '<', $info_path or throw("cannot open $info_path: $!");
        
        my $info_line = <$IN_info>; 
        print $OUT_INFO $info_line if $first_impute;      
      
        INFO_LINE:        
        while ($info_line = <$IN_info>) {
          my ($pos) = $info_line =~ /^\S+\s+\S+\s+(\d+)/; 
          next INFO_LINE if $pos < $bp_start->[$i];
          next INFO_LINE if $pos > $bp_end->[$i];
          print $OUT_INFO $info_line;
        }
        close $IN_info;
        $first_impute = 0; 
      }
      close $OUT_INFO;
   }
   
   ####
    $self->dbc->disconnect_when_inactive(0);

    $self->output_param('haps', $output_haps);
    $self->output_param('samples',  $output_sample);
    $self->output_param('impute', $output_impute) if $self->param_is_defined('impute');
    $self->output_param('impute_summary', $output_summary) if $self->param_is_defined('impute_summary');
    $self->output_param('impute_info', $output_info) if $self->param_is_defined('impute_info');
    

}
1;

