package ReseqTrack::Tools::RunPeakCall;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunProgram;

use Data::Dumper;
use File::Basename;
use base qw(ReseqTrack::Tools::RunProgram);

=pod

=head1 NAME

ReseqTrack::Tools::RunPeakCall::Fseq

=head1 SYNOPSIS

Base class for running peak callers. Child class of ReseqTrack::Tools::RunProgram

Implementations should be able to deal with a bed file, and optionally a bam file.
This class will convert bam to bed if required
If the BAM file has duplicates marked, the module will remove them prior to conversion if strip_duplicates is set to a value that evaluates to true

=head1 Example

my $fseq = ReseqTrack::Tools::RunPeakCall::Fseq(
					-input_files => '/path/to/file.bed'
					-working_dir => $output_dir,
					-options => \%options,
					-job_name => $name,
					-strip_duplicates => 1,
                      );
$fseq->run;
my $output_file_list = $fseq->output_files;

=cut

sub new {
	my ( $class, @args ) = @_;
	
	my $self = $class->SUPER::new(@args);

	my ( $bam_to_bed_path, $strip_duplicates, $samtools_path, $control_files )
	    = rearrange( [
	         qw( BAM_TO_BED_PATH STRIP_DUPLICATES
				SAMTOOLS_PATH CONTROL_FILES)
			], @args);
	
	$self->bam_to_bed_path($bam_to_bed_path || 'bamToBed');	
	$self->samtools_path($samtools_path || 'samtools');	
	$self->strip_duplicates($strip_duplicates);
	$self->control_files($control_files);
	
	if ($self->control_required && ! scalar(@{$self->control_files})){
		throw("Control files required but not specified");
	}
	if (scalar(@{$self->control_files}) && ! $self->can_use_control){
		thow("Control files specified but cannot be used");
	}

	return $self;
}

sub can_read_bam {
	throw("can_read_bam should be implemented in the child class");
}

sub can_use_control {
	throw("can_use_control should be implemented in the child class");
}
sub control_required {
	throw("control_required should be implemented in the child class");
}


sub convert_file {
	my ($self, $file_name) = @_;
	
	my $bam_to_bed_path = $self->bam_to_bed_path;
	my $samtools_path = $self->samtools_path;
	
	my $temp_dir = $self->get_bkilltemp_dir;
	my ($name,$path,$suffix) = fileparse($file_name);
	
	my $target_name = $temp_dir.'/'.$name;
	
	my $cmd;
	
	if ($file_name =~ m/\.bam$/ ){
		if ($self->can_read_bam){
			if ($self->strip_duplicates){
				$cmd = "$samtools_path view -F 1024 -F 4 -bh $file_name > $target_name";
			}
			else {
				$cmd = undef;
			}
		}
		else{
			if ($self->strip_duplicates){
				$cmd = "$samtools_path view -F 1024 -F 4 -buh $file_name | $bam_to_bed_path -i - > $target_name";
			}
			else {
				$cmd = "$bam_to_bed_path -i $file_name > $target_name";
			}			
		}
	}
	
	if ($cmd) {
		$self->execute_command_line($cmd);
		return $target_name;
	}
	else {
		return $file_name;
	}
}

sub get_input_files_cmd{
	my ($self) = @_;
	my @converted_files = map {$self->convert_file($_)} @{$self->input_files};
	return \@converted_files; 
}

sub get_control_files_cmd{
	my ($self) = @_;
	my @converted_files = map {$self->convert_file($_)} @{$self->control_files};
	return \@converted_files; 
}


=head2 control_files

  Arg [1]   : ReseqTrack::Tools::RunPeakCall
  Arg [2]   : string or arrayref of strings
  Function  : accessor method for control files. Confusingly, these are often called 'input' files (chromatin input, not enriched for ChIP target)
  Returntype: arrayref of strings
  Example   : $self->control_files('path/to/file');

=cut

sub control_files {
  my ( $self, $arg ) = @_;

  $self->{'control_files'} ||= {};
  if ($arg) {
    foreach my $file (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      $file =~ s{//}{/}g;
      $self->{'control_files'}->{$file} = 1;
    }
  }

  my @files = keys %{$self->{'control_files'}};
  return \@files;
}

sub bam_to_bed_path {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'bam_to_bed_path'} = $arg;
    }
    return $self->{'bam_to_bed_path'};
}

sub samtools_path {
    my ($self, $arg) = @_;
    if ($arg) {
        $self->{'samtools_path'} = $arg;
    }
    return $self->{'samtools_path'};
}

sub strip_duplicates {
    my ($self, $arg) = @_;
    if (defined $arg) {
        $self->{'strip_duplicates'} = $arg;
    }
    return $self->{'strip_duplicates'};
}

sub bed_file {
    my ($self, $arg) = @_;
    if (defined $arg) {
        $self->{'bed_file'} = $arg;
    }
    return $self->{'bed_file'};
}

1;