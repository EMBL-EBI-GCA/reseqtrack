package ReseqTrack::Tools::GATKTools::IndelReAligner;

use strict;
use warnings;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::AlignmentBase;
use ReseqTrack::Tools::Argument qw(rearrange);
use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::GATKTools);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ($intervals_file) = rearrange(
		[
			qw(
			  INTERVALS_FILE
			  )
		],
		@args
	);

	$self->intervals_file($intervals_file);
	$self->gatk_tool("Indelrealigner");
	$self->check_jar_file_exists();
	
	if ( !defined $self->options ) {
		#print $self->gatk_tool ." Using default options\n";
		$self->options(
			" -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing ");
	}

	return $self;
}

sub run {

	my $self = shift;

	my $cmd = $self->construct_run_cmd();

	#$self->execute_command_line( $cmd );

	return;

}

sub construct_run_cmd {
	my $self = shift;

	my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
	$cmd .= $self->gatk_path . "\/" . $self->jar_file;
	$cmd .= " -T " . $self->gatk_tool . " ";

	$cmd .= $self->options;
	$cmd .= "-targetIntervals " . $self->intervals_file . " ";

	my $knowns = $self->known;
	foreach my $k (@$knowns) {
		$cmd .= "-known $k ";
	}

	$cmd .= " -R " . $self->reference . " ";
	$cmd .= " -I " . $self->bam;

	my $realigned_bam = $self->bam . ".indel_realigned.bam";
	$cmd .= " -o " . $realigned_bam;
	print "\n\n$cmd\n";

	$self->output_files($realigned_bam);
	
	#$self->files_to_delete($self->intervals_file);
	#$self->files_to_delete($self->bam);

	return $cmd;

}

sub intervals_file {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{intervals_file} = $arg;
	}
	return $self->{intervals_file};
}

1;
