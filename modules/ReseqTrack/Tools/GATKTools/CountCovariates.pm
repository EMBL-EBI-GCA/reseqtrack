package ReseqTrack::Tools::GATKTools::CountCovariates;

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

	my ($known_sites) = rearrange(
		[
			qw(
			  KNOWN_SITES
			  )
		],
		@args
	);

	my $default_options = "";
	$default_options .= " -l INFO -L '1;2;3;4;5;6;7;8;9;10;11;12;";
	$default_options .= "13;14;15;16;17;18;19;20;21;22;X;Y;MT' ";
	$default_options .= "-cov ReadGroupCovariate -cov QualityScoreCovariate ";
	$default_options .= "-cov CycleCovariate -cov DinucCovariate ";

	
	$self->gatk_tool("CountCovariates");
	$self->known_sites($known_sites);
	$self->check_jar_file_exists();

	if ( !defined $self->options ) {
			#print $self->gatk_tool .": Using default options\n";
			$self->options($default_options);
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

	$cmd .= $self->options . " ";

	$cmd .= " -R " . $self->reference . " ";
	$cmd .= " -I " . $self->bam . " ";

	my $known_sites = $self->known_sites;

	foreach my $k_sites (@$known_sites) {
		$cmd .= " -knownSites $k_sites ";
	}

	my $covariates_file = $self->bam . "\.recal_data.csv ";
	$cmd .= " -recalFile $covariates_file";

	$self->output_files($covariates_file );

	print "\n\n$cmd\n";

	return $cmd;
}



1;
