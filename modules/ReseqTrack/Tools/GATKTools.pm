
=pod

=head1 NAME

ReseqTrack::Tools::GATKTools

=head1 SYNOPSIS

Object to create a bam file that in realigned around
known indel sites using GATK RealignerTargetCreator and
IndelRealigner

ISA = ReseqTrack::Tools::RunAlignment

If ENV's $SAMTOOLS and $GATK are set, then no need to 
pass to object. They will be assigned accordingly.

SAMTOOLS required to index bams and extract headers.

example

my $GATK = $IR->new(
		     -java_exe        =>"/usr/bin/java" ,
		     -jvm_args        =>"-Xmx4g",
 		     -GATK_PATH     =>$GATK/GenomeAnalysisTK/",
		     -working_dir     => $input{working_dir},
);



=cut


package ReseqTrack::Tools::GATKTools;

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;

use ReseqTrack::Tools::RunProgram;
use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::RunProgram);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $java_exe, $jvm_args, $options, $gatk_path ,
	     $reference) = rearrange(
		[
			qw(
			  JAVA_EXE
			  JVM_ARGS
			  OPTIONS
			  GATK_PATH
REFERENCE
			  )
		],
		@args
	);

	#defaults
	$self->java_exe("/usr/bin/java");
	$self->jvm_args("-Xmx4g");
	$self->reference($reference);
	$self->java_exe($java_exe);
	$self->jvm_args($jvm_args);
	$self->gatk_path($gatk_path);


	if ( ! defined ($self->samtools)){

	  if ( defined $ENV{SAMTOOLS}){
	    my $exe = $ENV{SAMTOOLS} . '/'. "samtools";
	    $exe =~ s|//|\/|g;
	    print "Setting path to samtools from ENV\n";
	    $self->samtools($exe) if ( -e $exe);

	    if (! -e $exe){
	      throw "Path to samtools not set. Required for indexing bams\n"; 
	    }

	  }
	}

	if (!defined $self->gatk_path){
	  if (defined $ENV{GATK}){
	    print "Setting GATK path from ENV\n";
	    	$self->gatk_path($ENV{GATK});
	  }
	  else{
	    throw "Path to GATK directory not set\n";
	  }
	}

	$self->bam(${$self->input_files}[0]);
	return $self;
}



sub check_bai_exist{
  my ($self,$bam) = @_;

  my $samtools = $self->samtools; 
  my $bamindex = $bam . "\.bai";

  if ( !-e $bamindex){
    print "$bamindex does not exist. Creating\n";
  }
  else{
    return;
  }

  my $cmd = "$samtools index $bam";

  `$cmd`;

  $self->files_to_delete($bamindex);

  print "Created $bamindex\n\n";
  return;

}


=head2 options

  Arg [1]   : ReseqTrack::Tools::CallBySamtools or CallByUmake or CallByGATK
  Arg [2]   : string, name of key e.g. "mpileup"
  Arg [3]   : string, optional, options to be used on command line e.g. "-m 100000000"
  Function  : accessor method for command line options
  Returntype: string, command line options
  Exceptions: n/a
  Example   : my $mpileup_options = $self->options{'mpileup'};

=cut

sub options {
    my ($self, $option_name, $option_value) = @_;

    if (! $self->{'options'}) {
        $self->{'options'} = {};
    }

    throw( "option_name not specified")
        if (! $option_name);

    if ($option_value) {
        $self->{'options'}->{$option_name} = $option_value;
    }

    return $self->{'options'}->{$option_name};
}





sub check_jar_file_exists {
	my ($self) = shift;

	my $jar_file = $self->gatk_path . "\/" . $self->jar_file;
	throw "Could not find:$jar_file" if ( !-e $jar_file );
	return;
}

sub java_exe {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{java_exe} = $arg;
	}
	return $self->{java_exe};
}

sub jvm_args {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{jvm_args} = $arg;
	}
	return $self->{jvm_args};
}


sub gatk_path {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{gatk_path} = $arg;
	}
	return $self->{gatk_path};

}

sub jar_file {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{jar_file} = $arg;
	}
	return $self->{jar_file};

}

sub gatk_tool {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{gatk_tool} = $arg;
	}
	return $self->{gatk_tool};

}


sub bam {
	my ( $self, $arg ) = @_;
	if ($arg) {
	  throw "Bam does not exist:$arg\n" if (! -e $arg); 
		$self->{bam} = $arg;
	}
	return $self->{bam};

}

sub known {
  my ( $self, $arg ) = @_;
  print "adding $arg\n" if ($arg);
  if ($arg) {
    if ( ref($arg) eq 'ARRAY' ) {
      foreach my $file (@$arg) {
	check_file_exists($file);
	push( @{ $self->{'known'} }, $file );
      }
    } else {
      check_file_exists($arg);
      push( @{ $self->{'known'} }, $arg );
    }
  }
	return $self->{known};
}

sub known_sites {
	my ( $self, $arg ) = @_;
	print "adding sites  $arg\n" if ($arg);
	
	if ($arg) {
		if ( !-e $arg ) {
			throw "Known sites file: $arg does not exist\n";
		}
		push( @{ $self->{known_sites} }, $arg );
	}

	return $self->{known_sites};

}

sub samtools {
  my ( $self, $arg ) = @_;
	
  if ($arg) {
    if ( !-e $arg ) {
      throw "Samtools exe: $arg does not exist\n";
    }
     $self->{samtools} = $arg;
  }

  return $self->{samtools};

}
sub reference {
  my ( $self, $arg ) = @_;
	
  if ($arg) {
    if ( !-e $arg ) {
      throw "Reference file: $arg does not exist\n";
    }
    $self->{reference} = $arg;
  }

  return $self->{reference};

}

1;
