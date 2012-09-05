package ReseqTrack::Tools::RunAlignment::GSNAP;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use Data::Dumper;
use File::Basename;

@ISA = qw(ReseqTrack::Tools::RunAlignment);

=pod

=head1 NAME

ReseqTrack::Tools::RunAlignment::GSNAP

=head1 SYNOPSIS

class for running GSNAP. Child class of ReseqTrack::Tools::RunAlignment

=head1 Example

my $gsnap = ReseqTrack::Tools::RunAlignment::GSNAP(
                      -input => '/path/to/file'
					  -map_dir => '/path/to/references',
                      -reference => 'name of reference',
                      -options => {'threads' => 4},
                      -read_group_fields => {'ID' => 1, 'LB' => 'my_lib'},
                      -first_read => 1000,
                      -last_read => 2000,
                      );
$gsnap->run;
my $output_file_list = $gsnap->output_files;

=cut

sub DEFAULT_OPTIONS { return {
		'batch_mode' => 5,
		'threads' => 1,
		'known_snps_file' => undef,
		'known_splice_sites_file' => undef,
		'novel_splice_sites' => 0,
		'split_output' => 1, # use GSNAPS split output. The good quality fraction (uniq, mult) will be merged into the output file
		'kmer_size' => undef,
		'orientation' => 'FR',
		'max_mismatches' => undef,
		'indel_penalty' => undef,
		'trim_mismatch_score' => undef,
		'trim_indel_score' => undef,
		'distant_splice_penalty' => undef,
		'local_splice_dist' => undef,
		'max_multihits' => undef,
		'nofails'	=> undef,
		'quiet_if_excessive' => undef,
		'map_dir' => undef,
	};
}

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	#setting defaults
	if (!$self->program) {
		if ($ENV{GSNAP}) {
			$self->program($ENV{GSNAP} . '/gsnap');
		}
		else {
			$self->program('gsnap');
		}
	}

	return $self;
}

sub run_alignment {
    my ($self) = @_;
    
	my $output_file = $self->working_dir() . '/'. $self->job_name;
	$output_file =~ s{//}{/};
	
	
	if ($self->fragment_file){
		$self->_do_alignment($output_file.'_se',$self->fragment_file);
	}
	if ( $self->mate1_file && $self->mate2_file ) {
		$self->_do_alignment($output_file.'_se',$self->mate1_file, $self->mate2_file);
	 }
}

#
sub _do_alignment {
	my ($self,$output_file,@input_files) = @_;

	#GSNAP doesn't output BAM, but we can pipe the output through samtools
	my $do_bam_conversion = 0;
	my $output_option =  $self->output_format;
	if (lc($output_option) eq 'bam') {
		$output_option = 'sam';
		$do_bam_conversion = 1;
	}

	my @cmd_words = ("bash -c '");

	push(@cmd_words, $self->program);
	
	#reference information
	push(@cmd_words,'-D',$self->options('map_dir')) if ($self->options('map_dir'));
	push(@cmd_words,'-d',$self->reference);
	
	#parallel running
	push(@cmd_words,'-q',$self->batch_index.'/'.$self->batch_job_count) if ($self->batch_job_count);
	
	#output options - format and read group information
	push(@cmd_words, '-A', $output_option) if ($output_option);
	
	my $read_group_fields = $self->read_group_fields;
	if ($read_group_fields){
		my %fields_to_args = (
			'ID' => 'read-group-id',
			'SM' => 'read-group-name',
			'LB' => 'read-group-library',
			'PL' => 'read-group-platform',			
		);
		while (my ($field, $value) = each %$read_group_fields) {
			next unless $fields_to_args{$field};
			push (@cmd_words, '--'.$fields_to_args{$field}.'='.$self->escape_spaces($value));
		}
	}
		
	#operations - threads and memory modes	
	push(@cmd_words, '-B', $self->options('batch_mode'))
		if ($self->options('batch_mode'));
	push(@cmd_words, '-t', $self->options('threads'))
		if ($self->options('threads'));
	
	#splicing and snps 
	push(@cmd_words, '-s', $self->options('known_splice_sites_file'))
		if ($self->options('known_splice_sites_file'));
	push(@cmd_words, '-v', $self->options('known_snps_file'))
		if ($self->options('known_snps_file'));
	push(@cmd_words, '-N', $self->options('novel_splice_sites'))
		if ($self->options('novel_splice_sites'));

		
	# inputs
	if ($input_files[0] =~ /.gz$/){
		push(@cmd_words, '--gunzip');
	}
	
	push(@cmd_words, '-k', $self->options('kmer_size')) if ($self->options('kmer_size'));
	push(@cmd_words, '-o', $self->options('orientation')) if (defined $self->options('orientation'));
	push(@cmd_words, '-m', $self->options('max_mismatches')) if (defined $self->options('max_mismatches'));
	push(@cmd_words, '-i', $self->options('indel_penalty')) if (defined $self->options('indel_penalty'));
	push(@cmd_words, '--trim-mismatch-score='.$self->options('trim_mismatch_score')) if (defined $self->options('trim_mismatch_score'));
	push(@cmd_words, '--trim-indel-score='.$self->options('trim_indel_score')) if (defined $self->options('trim_indel_score'));
	push(@cmd_words, '-E', $self->options('distant_splice_penalty')) if (defined $self->options('distant_splice_penalty'));
	push(@cmd_words, '-w', $self->options('local_splice_dist')) if (defined $self->options('local_splice_dist'));
	push(@cmd_words, '-n', $self->options('max_multihits')) if (defined $self->options('max_multihits'));
	push(@cmd_words, '--nofails') if (defined $self->options('nofails'));
	push(@cmd_words, '--quiet-if-excessive') if ($self->options('quiet_if_excessive'));

	#output
	my $output_base;
	if ($self->options('split_output')){
		$output_base = $self->working_dir().'/'.$self->job_name;
		push(@cmd_words, '--split-output='.$output_base);
	} 
	
	push(@cmd_words, @input_files);
	
	if (! $self->options('split_output')) {
		# single output file
		if ($do_bam_conversion){
			push(@cmd_words, '|', $self->samtools, 'view -bhS -');
		}
	
		push(@cmd_words, '>', $output_file);
		$self->output_files($output_file); 
	}
	
	#closing quote for bash -c
	push(@cmd_words, "'");

	my $cmd = join(' ', @cmd_words);

    $self->execute_command_line($cmd);

	if ($self->options('split_output')){
		my @output_files;
		
		
		my @sam_files = glob($output_base.'*');
				
		for my $sam_file (@sam_files) {
			if (empty_sam_file($sam_file)){ # remove empty files
				unlink($sam_file);
			}
			elsif ($do_bam_conversion && $sam_file !~ m/\.bam$/){
				my $bam_file = $sam_file .'.bam';
				my $convert_command = $self->samtools." view -bhS $sam_file > $bam_file";
				# can fail here if file only contains header
				$self->execute_command_line($convert_command);
				unlink($sam_file);
				
				push @output_files, $bam_file;
			}
			else{
				push @output_files, $sam_file;
			}
		}
		
		$self->output_files(\@output_files);
	}
}

sub empty_sam_file{
	my ($sam_file) = @_;
	return 0 if ($sam_file =~ m/\.bam$/);
	return 1 if (-z $sam_file); # file contains zero bytes
	
	open my $SAM, '<', $sam_file or throw("could not read $sam_file: $!");
	while (<$SAM>){		
		if (! m/^@/){
			return 0; # file contains non-header line
		}
	}
	close $SAM;
	#file only contains header
	return 1;
	
}

sub batch_index {
	my ($self) = @_;
	return $self->first_read;
}

sub batch_job_count {
	my ($self) = @_;
	return $self->last_read;
}

sub escape_spaces {
	my ($self,$string) = @_;
	$string =~ s/ /\\ /g;
	return $string;
}

1;