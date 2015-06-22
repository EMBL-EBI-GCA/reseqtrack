package ReseqTrack::Hive::Process::MergeByPaste;

use strict;
use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists get_lines_from_file);
use File::Basename qw(basename);
use ReseqTrack::Tools::Concatenate_by_paste;

sub param_defaults {
  return {
	program => undef,
	site_vcf => undef,
    save_files_from_deletion => undef,  
  };
}

sub run {
    my $self = shift @_;

    $self->param_required('files');
    my $files =  $self->param_as_array('files');
    my $program = $self->param('program');
    my $save_files_from_deletion =  $self->param('save_files_from_deletion');
    
    #my @sorted_files = sort{$a cmp $b} @$files; ### No need to sort here as the $self->input_files random the order again
       
    my $site_vcf = $self->param('site_vcf');
    check_file_exists($site_vcf);
	my $site_file_lines = get_lines_from_file($site_vcf);
	my $sites = $self->output_dir . "/" . basename($site_vcf) . ".$$.sites";  ### the file is written to #vcf_base_name# level, not to the nested #sample# level
	open (OUT, ">", $sites) || throw("Cannot open output file to write sites $sites");
	foreach my $line (@$site_file_lines) {
		next if $line =~ /^\#\#/;
		my @data = split (/\t/, $line);
		print OUT "$data[0]\t$data[1]\n";
	}
	close(OUT);
	
	unshift @$files, $sites;
	
	my $object = ReseqTrack::Tools::Concatenate_by_paste->new (
		-program					=> $program,
		-input_files				=> $files,
		-working_dir				=> $self->output_dir,
		-job_name					=> $self->job_name,
		-save_files_from_deletion	=> $save_files_from_deletion,
	);

	$self->run_program($object);
    $self->output_param('matrix', $object->output_files->[0]);
	unlink($sites);
}

1;
