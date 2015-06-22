package ReseqTrack::Hive::Process::GtVcf2SiteVcf;

use strict;
use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_executable);
use File::Basename qw(basename);
use IPC::System::Simple qw(system);

sub param_defaults {
  return {  
    tabix => undef,
  };
}

sub run {
    my $self = shift @_;
    my $gt_vcf = $self->param_required('gt_vcf');
	my $output_dir = $self->output_dir;

    check_file_exists($gt_vcf);   
        
    my $output_file = $output_dir . "/" . $self->job_name . ".site.vcf"; 
    
    my $command;
	if ($gt_vcf =~ /.vcf.gz$/) {
		$command = "zcat $gt_vcf | cut -f1-8 > $output_file";
	}
	else {
		$command = "cut -f1-8 $gt_vcf > $output_file";
	}	
    system($command);
    my $exit = $?>>8;
    throw("cut command $command died $!") if ($exit >1);
    
    $self->output_param('site_vcf', $output_file);
}

1;
