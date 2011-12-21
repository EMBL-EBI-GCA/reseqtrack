package ReseqTrack::Tools::GATKTools::TableRecalibration;

use strict;
use warnings;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunAlignment;
use Data::Dumper;

use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::GATKTools);

sub new {
	my ( $class, @args ) = @_;
	
	
	my $self = $class->SUPER::new(@args);

	my ( $known_sites ) = rearrange(
		[
			qw(
			  KNOWN_SITES
			  )
		],
		@args
	);

	$self->jar_file ("GenomeAnalysisTK.jar");
	$self->check_jar_file_exists();

	my $default_covariate_options = "";
	$default_covariate_options .= " -l INFO -L '1;2;3;4;5;6;7;8;9;10;11;12;";
	$default_covariate_options .= "13;14;15;16;17;18;19;20;21;22;X;Y;MT' ";
	$default_covariate_options .= "-cov ReadGroupCovariate -cov QualityScoreCovariate ";
	$default_covariate_options .= "-cov CycleCovariate -cov DinucCovariate ";

	$self->options ('table_recalibration_options', " -l INFO -compress 0 --disable_bam_indexing ");
	$self->options ('covariate_options'  ,$default_covariate_options);



    	$self->options ('knowns_sites'  ,$known_sites);
    
	if ( ! defined $self->{options}->{'knowns_sites'}){
	  throw ("No known sites defined");
	}

    
	$self->construct_make_covariates_file_cmd();
	$self->construct_table_recalibration_cmd();
	
	return $self;
}


sub run {

	my $self = shift;

	$self->create_covariates_file();

	$self->create_recalibrated_bam_file();

	return;

}



sub create_recalibrated_bam_file {
  my ( $self) = shift;
  
  $self->check_bai_exist($self->bam);
  $self->execute_command_line ($self->make_recalibrated_bam_cmd);

  $self->bam($self->recalibrated_bam_file_name);

  return;
}



sub create_covariates_file {
  my ( $self) = shift;
  
  $self->check_bai_exist($self->bam);
  $self->execute_command_line ($self->make_covariates_file_cmd);

  return;
}


sub construct_table_recalibration_cmd {

	my $self = shift;

	my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
	$cmd .= $self->gatk_path . "\/" . $self->jar_file;
	$cmd .= " -T TableRecalibration ";

	$cmd .=  $self->{options}->{'table_recalibration_options'};

	$cmd .= " -R " . $self->reference . " ";
	$cmd .= " -recalFile " . $self->covariates_file . " ";
	$cmd .= " -I " . $self->bam;

	my $recal_bam = $self->bam . "\.recal\.bam";
	$cmd .= " -o $recal_bam";

	$self->make_recalibrated_bam_cmd($cmd);

	$self->recalibrated_bam_file_name($recal_bam);

	print "\n\n$cmd\n\n";

	return $cmd;
}


sub construct_make_covariates_file_cmd {

	my $self = shift;

	my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
	$cmd .= $self->gatk_path . "\/" . $self->jar_file;
	$cmd .= " -T CountCovariates ";



	$cmd .= " -R " . $self->reference . " ";
	$cmd .= " -I " . $self->bam . " ";

	
	$cmd .=  " -knownSites ". $self->{options}->{'knowns_sites'} ." ";
		  

	my $covariates_file = $self->bam . "\.recal_data.csv";
	$cmd .= " -recalFile $covariates_file";

	$self->make_covariates_file_cmd ($cmd);

	$self->covariates_file($covariates_file );

	$self->files_to_delete ($covariates_file);

	print "\n\n$cmd\n\n\n";

	return;

}


sub make_recalibrated_bam_cmd {

  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{make_recalibrated_bam_cmd} = $arg;
  }

  return $self->{make_recalibrated_bam_cmd};
}


sub make_covariates_file_cmd {

  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{make_covariates_file_cmd} = $arg;
  }

  return $self->{make_covariates_file_cmd};
}



sub create_table_recalibration_cmd {

  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{create_table_recalibration_cmd} = $arg;
  }

  return $self->{create_table_recalibration_cmd};
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


sub recalibrated_bam_file_name {
  my ( $self, $arg ) = @_;

  if ($arg) {
    if ( !-e $arg ) {
      #warning "Covariates file: $arg does not exist\n";
    }
    $self->{recalibrated_bam_file_name} = $arg;
  }

  return $self->{recalibrated_bam_file_name};

}

1;
