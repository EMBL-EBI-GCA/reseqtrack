package ReseqTrack::Tools::RunPeakCall::Macs;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Tools::GeneralUtils qw(execute_system_command execute_pipe_system_command);

use Cwd;
use File::Copy;
use base qw(ReseqTrack::Tools::RunPeakCall);

=pod

=head1 NAME

ReseqTrack::Tools::RunPeakCall::Macs

=head1 SYNOPSIS

class for running MACS 1.4.2. Child class of ReseqTrack::Tools::RunPeakCall
http://fureylab.web.unc.edu/software/fseq/

Set fragment_size to 0 for DNAse data

Takes either a bed file or a bam file.
BAM files will be converted to bed with bamToBed
If the BAM file has duplicates marked, the module will remove them prior to conversion if strip_duplicates is set to a value that evaluates to true.

=head1 Example

my $macs = ReseqTrack::Tools::RunPeakCall::Macs(
					-input_files => '/path/to/file.bed'
					-working_dir => $output_dir,
					-options => \%options,
					-job_name => $name,
                      );
$macs->run;
my $output_file_list = $macs->output_files;

=cut
sub DEFAULT_OPTIONS { return {
	genome_size => 'hs', #  -g/--gsize #(hs) genome size
	tag_size => undef, #	-s/--tsize#	The size of sequencing tags. If you DON'T specify it, MACS will try to use the first 10 sequences from your input treatment file to determine the tag size. Specifying it will override the automatic determined tag size.
	bandwidth => undef, #	--bw#	The band width which is used to scan the genome for model building. You can set this parameter as the sonication fragment size expected from wet experiment. The previous side effect on the peak detection process has been removed. So this parameter only affects the model building.
	pvalue => undef, #-p/--pvalue #	The pvalue cutoff. Default is 1e-5.
	mfold => undef, #-m/--mfold #	This parameter is used to select the regions within MFOLD range of high-confidence enrichment ratio against background to build model. The regions must be lower than upper limit, and higher than the lower limit of fold enrichment. DEFAULT:10,30 means using all regions not too low (>10) and not too high (<30) to build paired-peaks model. If MACS can not find more than 100 regions to build model, it will use the --shiftsize parameter to continue the peak detection.
	nolambda => 0, #--nolambda #	With this flag on, MACS will use the background lambda as local lambda. This means MACS will not consider the local bias at peak candidate regions.
	slocal => undef,
	llocal => undef, #--slocal, --llocal#	These two parameters control which two levels of regions will be checked around the peak regions to calculate the maximum lambda as local lambda. By default, MACS considers 1000bp for small local region(--slocal), and 10000bps for large local region(--llocal) which captures the bias from a long range effect like an open chromatin domain. You can tweak these according to your project. Remember that if the region is set too small, a sharp spike in the input data may kill the significant peak.
	on_auto => 0, #--on-auto # 	Whether turn on the auto paired-peak model process. If set, when MACS failed to build paired model, it will use the nomodel settings, the '--shiftsize' parameter to shift and extend each tags. If not set, MACS will be terminated if paried-peak model is failed.
	nomodel => 0, #--nomodel #	While on, MACS will bypass building the shifting model.
	shiftsize => undef, #--shiftsize #	While '--nomodel' is set, MACS uses this parameter to shift tags to their midpoint. For example, if the size of binding region for your transcription factor is 200 bp, and you want to bypass the model building by MACS, this parameter can be set as 100. This option is only valid when --nomodel is set or when MACS fails to build paired-peak model.
	keep_dup => undef, # --keep-dup #	It controls the MACS behavior towards duplicate tags at the exact same location -- the same coordination and the same strand. The default 'auto' option makes MACS calculate the maximum tags at the exact same location based on binomal distribution using 1e-5 as pvalue cutoff; and the 'all' option keeps every tags. If an integer is given, at most this number of tags will be kept at the same location. Default: 1.
	to_large => undef, #--to-large #	When not set, scale the larger dataset down to the smaller dataset; when set, the smaller dataset will be scaled towards the larger dataset.
	wiggle => 0, #-w/--wig #	If this flag is on, MACS will store the fragment pileup in wiggle format for every chromosome. The gzipped wiggle files will be stored in subdirectories named NAME+'_MACS_wiggle/treat' for treatment data and NAME+'_MACS_wiggle/control' for control data. --single-profile option can be combined to generate a single wig file for the whole genome.
	bedgraph => 0, #-B/--bdg # 	If this flag is on, MACS will store the fragment pileup in bedGraph format for every chromosome. The bedGraph file is in general much smaller than wiggle file. However, The process will take a little bit longer than -w option, since theoratically 1bp resolution data will be saved. The bedGraph files will be gzipped and stored in subdirectories named NAME+'_MACS_bedGraph/treat' for treatment and NAME+'_MACS_bedGraph/control' for control data. --single-profile option can be combined to generate a single bedGraph file for the whole genome.
	space => undef, #--space=SPACE #	By default, the resoluation for saving wiggle files is 10 bps,i.e., MACS will save the raw tag count every 10 bps. You can change it along with '--wig' option.
	call_subpeaks => 0, #--call-subpeaks #	If set, MACS will invoke Mali Salmon's PeakSplitter software through system call. If PeakSplitter can't be found, an instruction will be shown for downloading and installing the PeakSplitter package. The PeakSplitter can refine the MACS peaks and split the wide peaks into smaller subpeaks. For more information, please check the following URL:
	verbose => undef,
	};
}

sub new {
	my ( $class, @args ) = @_;
	
	my $self = $class->SUPER::new(@args);

	#setting defaults
	if (!$self->program) {
	   if ($ENV{fseq}) {
	     $self->program($ENV{macs} . '/macs14');
	   }
	   else {
	     $self->program('macs14');
	   }
	 }

	return $self;
}

sub can_read_bam {
	return 1;
}

sub can_use_control {
	return 1;
}
sub control_required {
	return 0;
}

sub run_program {
	my ($self) = @_;
	
	my @cmd_args;
	my $job_name = $self->job_name;
	my $temp_dir = $self->get_temp_dir();
	
	push @cmd_args, $self->program;

	for my $input_file (@{$self->input_files}){
		push @cmd_args, '-t', $input_file;
	}

	push @cmd_args, '-n', $job_name;

	push @cmd_args, '--gsize', $self->options('genome_size') if ($self->options('genome_size')); 
	
	
	if (@{$self->control_files}){
	  push @cmd_args, '-c <(samtools merge - ';
	  for my $control_file (@{$self->control_files}){
  		push @cmd_args, $control_file;
  	}
	  push @cmd_args, ')';
	}
	#for my $control_file (@{$self->control_files}){
	#	push @cmd_args, '-c', $control_file;
	#}
	
	push @cmd_args, '--tsize', $self->options('tag_size') if ($self->options('tag_size')); 
	push @cmd_args, '--bw='.$self->options('bandwidth') if ($self->options('bandwidth')); 
	push @cmd_args, '--pvalue', $self->options('pvalue') if ($self->options('pvalue')); 
	push @cmd_args, '--mfold', $self->options('mfold') if ($self->options('mfold')); 
	push @cmd_args, '--nolambda' if ($self->options('nolambda')); 
	push @cmd_args, '--slocal='.$self->options('slocal') if ($self->options('slocal')); 
	push @cmd_args, '--llocal='.$self->options('llocal') if ($self->options('llocal')); 
	push @cmd_args, '--on-auto' if ($self->options('on_auto')); 
	push @cmd_args, '--nomodel' if ($self->options('nomodel'));
	push @cmd_args, '--shiftsize='.$self->options('shiftsize') if (defined $self->options('shiftsize'));
	push @cmd_args, '--keep-dup='.$self->options('keep_dup') if ($self->options('keep_dup')); 
	push @cmd_args, '--to-large' if ($self->options('to_large'));
	push @cmd_args, '--wig' if ($self->options('wiggle'));
	push @cmd_args, '--bdg' if ($self->options('bedgraph'));
	push @cmd_args, '--space='.$self->options('space') if ($self->options('space'));
	push @cmd_args, '--single-profile' if ($self->options('wiggle') || $self->options('bedgraph'));
	push @cmd_args, '--call-subpeaks' if ($self->options('call_subpeaks'));
	push @cmd_args, '--verbose='.$self->options('verbose') if (defined $self->options('verbose')); 
	
	my $dir = getcwd;
	# all output will go into a temp dir and be deleted by default. rescue the bits we need by moving them to the final output dir
	
	print "Changing dir to $temp_dir$/" if ($self->echo_cmd_line);
	chdir($temp_dir);
	
   	my $cmd = "bash -c '".join(' ', @cmd_args)."'";
	#print "Executing: $cmd";
	#execute_system_command($cmd); # this failed silently if run with $self->execute_command_line
	$self->execute_command_line(join(' ', @cmd_args));
	chdir($dir);
	print "Returning dir to $dir$/" if ($self->echo_cmd_line);# wrapping the cmd like this allows the <(..) file conversion to work
	
	my $root_temp_output = $temp_dir.'/'.$job_name;
	my $root_final_output = $self->working_dir. '/'.$job_name;
	
	my $bed_file_suffix  ='_peaks.bed';
	my @output_suffixes = ($bed_file_suffix,'_summits.bed', '_model.r', '_peaks.xls', '_negative_peaks.xls');

	for my $suffix (@output_suffixes){
		my $src = $root_temp_output.$suffix;
		my $dest = $root_final_output.$suffix;
		
		next unless (-e $src);
		
		# MACS uses the xls suffix for text output
		if ($suffix =~ m/xls$/){
			$dest =~ s/xls$/txt/;
		}
		
		print "Moving $src to $dest$/" if ($self->echo_cmd_line);
		
		move($src , $dest) or throw( "Failed to move $src to $dest: $!");
		$self->output_files($dest);
	}
	
	
	$self->bed_file($root_final_output.$bed_file_suffix);
	
	if ($self->options('wiggle')){
		my $wiggle_dir = $temp_dir.'_MACS_wiggle';
		my $dest = $root_final_output.'_treat_afterfitting.wig.gz'; 
		my $src = $wiggle_dir.'/'.$job_name.'_treat_afterfitting.wig.gz';
		move($src,$dest) or throw( "Failed to move $src to $dest ");
		$self->output_files($dest);
	}
	if ($self->options('bedgraph')){
		my $bdg_dir = $temp_dir.'_MACS_wiggle';
		my $dest = $root_final_output.'_treat_afterfitting.bdg.gz'; 
		my $src = $bdg_dir.'/'.$job_name.'_treat_afterfitting.bdg.gz';
		move($src,$dest) or throw( "Failed to move $src to $dest ");
		$self->output_files($dest);
	}

}

1;