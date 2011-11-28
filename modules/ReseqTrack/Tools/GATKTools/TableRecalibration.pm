package ReseqTrack::Tools::GATKTools::TableRecalibration;

use strict;
use warnings;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::AlignmentBase;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::GATKTools);

sub new {
	my ( $class, @args ) = @_;
	
	
	my $self = $class->SUPER::new(@args);

	my ( $covariates_file ) = rearrange(
		[
			qw(
			  COVARIATES_FILE
			  )
		],
		@args
	);

	$self->gatk_tool("TableRecalibration");
	$self->check_jar_file_exists();
    $self->covariates_file($covariates_file);
    
    
    my $default_options = " -l INFO -compress 0 --disable_bam_indexing ";
    
    if ( !defined $self->options ) {
		print  $self->gatk_tool .":Using default options\n";
		$self->options($default_options);
	}
	return $self;
}

sub run {

	my $self = shift;

	my $cmd = $self->construct_run_cmd();

	#$self->execute_command_line($cmd);

	return;

}

sub construct_run_cmd {

	my $self = shift;

	my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
	$cmd .= $self->gatk_path . "\/" . $self->jar_file;
	$cmd .= " -T " . $self->gatk_tool . " ";

	$cmd .= $self->options;
	$cmd .= " -R " . $self->reference . " ";
	$cmd .= " -recalFile " . $self->covariates_file . " ";
	$cmd .= " -I " . $self->bam;

	my $recal_bam = $self->bam . "\.recal\.bam";
	$cmd .= " -o $recal_bam";


	print "\n\n$cmd";

	return $cmd;
}

sub covariates_file {
	my ( $self, $arg ) = @_;

	if ($arg) {
		if ( !-e $arg ) {
			#warning "Covariates file: $arg does not exist\n";
		}
		$self->{covariates_file} = $arg;
	}

	return $self->{covariates_file};

}

=head

sub bam_recalibrate_quality_scores {
 
  $self->check_bai_present ($self->bam);


 

  my $recal_bam = $self->bam . "\.recal\.bam";
  $cmd .= " -o $recal_bam";

	
  $self->bam ( $recal_bam);
	
  print "\n\n$cmd\n";
  print "\nTableRecalibration Output = $recal_bam\n";
  eval{
    `$cmd`;
  };

  throw "TableRecalibration failed\n$@\n" if ($@);

  return ($recal_bam);
}

=cut

1;
