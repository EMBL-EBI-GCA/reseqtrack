package ReseqTrack::Tools::GATKTools::RealignerTargetCreator;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::RunProgram;

use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::GATKTools);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

=head	

	  my (   )
    = rearrange(
		[
		 qw(
			
			  )
		],
		@args
	       );

=cut	

	$self->gatk_tool("RealignerTargetCreator");
	$self->check_jar_file_exists();
	return $self;
}

sub run {
	my $self = shift;
	my $cmd = $self->construct_run_cmd();
	#$self->execute_command_line( $cmd );
	
	foreach (@{$self->output_files}){
		#throw ($self->gatk_tool ." :Does not exist:$_") if ( !-e $_);		
	}	
	return;
}

sub construct_run_cmd {
	
	my $self = shift;

	my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
	$cmd .= $self->gatk_path . "\/" . $self->jar_file;
	$cmd .= " -T " . $self->gatk_tool . " ";

	if ( defined $self->options ) {
		$cmd .= $self->options . " ";
	}
	else {
		my $knowns = $self->known;
		if ( !$knowns ) {
			throw("Must have known sites to realign around.None given");
		}
		foreach my $k (@$knowns) {
			$cmd .= "-known $k ";
		}
	}

	my $interval_file = $self->bam . ".interval_list ";
	
	$self->output_files ($interval_file);
	$cmd .= " -o $interval_file ";
	$cmd .= " -R " . $self->reference . " ";
	$cmd .= " -I " . $self->bam;

	print "\n$cmd\n";
	
	return $cmd;
}



1;


