package ReseqTrack::Hive::Process::VcfToBed;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw warning);
use File::Basename qw(fileparse);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_executable check_file_exists);


=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    $self->param_required('vcf');
    my $bgzip = $self->param_is_defined('bgzip') ? $self->param('bgzip') : 'bgzip';
    my $max_variants = $self->param_is_defined('max_variants') ? $self->param('max_variants') : 5000 ;
    my $overlap_variants = $self->param_is_defined('overlap_variants') ? $self->param('overlap_variants') : 0 ;
        
    my $vcfs = $self->param_as_array('vcf');

  
    throw("expecting single vcf") if(scalar @{$vcfs} > 1);
    
    foreach my $vcf_path (grep {defined $_} @$vcfs) {
      check_file_exists($vcf_path);
    }

    my $output_dir = $self->output_dir;
    my $job_name = $self->job_name;
    
    check_directory_exists($output_dir);
    check_executable($bgzip);
    
    my $chrom_start;
    my $chrom_end;
    
    my $output_file;
    my $SQ_start;
    my $SQ_end;
    my $vcf_name;
    
    
    my $first_vcf = 1;
    $self->dbc->disconnect_when_inactive(1);
    VCF:
    foreach my $i (0..$#{$vcfs}) {
      my $vcf_path = $vcfs->[$i];
      
      my $filename = fileparse($vcf_path, qw( .vcf .vcf.gz ));
      $vcf_name = $filename;
              
      my $output_bed = "$output_dir/$filename.bed";
      
      
               
      push @{$output_file}, $output_bed;
      open my $OUT, "> $output_bed";
      
      next VCF if ! defined $vcf_path;
      
      my $IN;
      if ($vcf_path =~ /\.b?gz(?:ip)?$/){
        open $IN, "$bgzip -cd $vcf_path |" or throw("cannot open $vcf_path: $!");
        
      }
      else {
        open $IN, '<', $vcf_path or throw("cannot open $vcf_path: $!");
      }
           
      my $start_chrom;
      my $start_pos;
      my $end_pos;
      my $chrom;
      my $pos;
      my $count = 0;
      my $overlap_pos;
      
      LINE:
      while (my $line = <$IN>) {
        if ($line =~ /^#/) {
            next LINE;
        }
        
        my @line_array = split("\t",$line);
           
        $count++;
        
        if($count == 1) {
          $start_chrom = $line_array[0];
          $start_pos = $line_array[1];
          $end_pos = $line_array[1];
          
          push @{$SQ_start}, $start_chrom;
             
        }
             
        $chrom = $line_array[0];
        $pos = $line_array[1];
        
        if($overlap_variants == $max_variants-$count+1) {
          $overlap_pos = $pos;
        }
        
        if($chrom eq $start_chrom) {
          if($count > $max_variants) {
            next if($pos == $end_pos);
          
            print $OUT "$chrom\t$start_pos\t$end_pos\t.\t1\t+\n";
            $count = $overlap_variants+1;
            
            $start_chrom = $chrom;
            $start_pos = $overlap_pos;
            $end_pos = $pos;
          }         
          else {
            $end_pos = $pos;
          }
        }
        else {
          print $OUT "$start_chrom\t$start_pos\t$end_pos\t.\t1\t+\n";
          $count = 1;
          
          $start_chrom = $chrom;
          $start_pos = $pos;
          $end_pos = $pos;
        }        
      }
      close $IN;
      
      print $OUT "$chrom\t$start_pos\t$pos\t.\t1\t+\n";
      
      push @{$SQ_end}, $chrom;
      
      $first_vcf = 0;
      
      close $OUT;
      
      $chrom_start= $$SQ_start[0];
      $chrom_end= $$SQ_end[0];
      
    }
    

    $self->dbc->disconnect_when_inactive(0);

    $self->output_param('bed', $$output_file[0]);
    $self->output_param('CHROM_start' , $chrom_start);
    $self->output_param('CHROM_end' , $chrom_end);
   # $self->output_param('vcf_name', $vcf_name);

}


1;

