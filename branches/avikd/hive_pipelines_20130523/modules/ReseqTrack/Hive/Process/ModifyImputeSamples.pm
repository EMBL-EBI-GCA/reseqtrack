package ReseqTrack::Hive::Process::ModifyImputeSamples;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use File::Basename qw(fileparse);
use ReseqTrack::Tools::FileSystemUtils qw(check_directory_exists check_executable check_file_exists);

=head2 run

    Description : Implements run() interface method of Bio::EnsEMBL::Hive::Process that is used to perform the main bulk of the job (minus input and output).

=cut

sub run {
    my $self = shift @_;
    
    $self->param_required('samples');
    $self->param_required('samples_gender_info');
    
    my $samples_gender_info = $self->param('samples_gender_info');
    
    check_file_exists($samples_gender_info);
    
    my %samples_gender_hash;
    
    open(my $INFO, '<', $samples_gender_info)or throw("cannot open $samples_gender_info: $!");
    
    while(my $info_line=<$INFO>){
    
      my($sample_name,$sex) = split('\s+', $info_line);
      $sex = 0 if(!$sex);
      $samples_gender_hash{$sample_name} = $sex;
    }
    close($INFO);
    
    my $samples = $self->file_param_to_flat_array('samples');
    
    throw("expecting single vcf") if(scalar @{$samples} > 1);
    
    foreach my $samples_path (grep {defined $_} @$samples) {
      check_file_exists($samples_path);
    }
    
    my $output_dir = $self->output_dir;
    my $job_name = $self->job_name;
    
    check_directory_exists($output_dir);
    
    my $output_file;
    
     SAMPLES:
    foreach my $i (0..$#{$samples}) {
    
      my $samples_path = $samples->[$i];
      my $filename = fileparse($samples_path, qw( .samples ));
      
      my $output_samples = "$output_dir/$filename.mod.samples";
      
      push @{$output_file}, $output_samples;
      
      open my $OUT, "> $output_samples";
      
      next SAMPLES if ! defined $samples_path;
      
      my $IN;
      open $IN, '<', $samples_path or throw("cannot open $samples_path: $!");
      
      my $in_line = <$IN>;
      chomp($in_line);
      print $OUT $in_line," sex\n";
      
      $in_line = <$IN>;
      chomp($in_line);
      print $OUT $in_line," D\n";
      
      while( $in_line = <$IN> ){
        chomp($in_line);
        
        if($in_line=~ /(\S+)\s+/){
          my $sample_name = $1 ; 
        
          if(exists($samples_gender_hash{$sample_name})) {
        
            if($samples_gender_hash{$sample_name} eq "M" ){
              print $OUT $in_line," 1\n";
            }
            elsif($samples_gender_hash{$sample_name} eq "F"){
              print $OUT $in_line," 2\n";
            }
            elsif($samples_gender_hash{$sample_name} eq "0"){
              print $OUT $in_line," 0\n";
            }
            else
            {
              throw("got gender: $samples_gender_hash{$sample_name} , expecting M / F");
            }
          }
          else
          {
            print $OUT $in_line," 0\n";
          }
        }
        else
        {
          throw("no sample name found");
        }
      }
      close $IN;
      close $OUT;  
    }
    
    $self->output_param('samples' , $output_file);
      
  }

1;