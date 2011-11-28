package ReseqTrack::Tools::GATKTools;

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::RunAlignment;
use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::RunAlignment);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $java_exe, $jvm_args, $options, $gatk_path, $jar_file ) = rearrange(
		[
			qw(
			  JAVA_EXE
			  JVM_ARGS
			  OPTIONS
			  GATK_PATH
			  JAR_FILE
			  )
		],
		@args
	);

	#defaults
	$self->java_exe("/usr/bin/java");
	$self->jvm_args("-Xmx4g");
	#$self->samtools("/home/smithre/Work/bin/samtools/samtools");
	$self->gatk_path("/home/smithre/Work/bin/GenomeAnalysisTK-1.2-61-g86871bd");

	$self->options($options);
	$self->java_exe($java_exe);
	$self->jvm_args($jvm_args);
	$self->gatk_path($gatk_path);
	$self->jar_file($jar_file);

	print Dumper $self;
	exit;

	return $self;
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

sub options {
	my ( $self, $arg ) = @_;
	if ($arg) {
		$self->{options} = $arg;
	}
	return $self->{options};

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

sub known {
	my ( $self, $arg ) = @_;
	print "adding $arg\n" if ($arg);
	if ($arg) {
		if ( !-e $arg ) {
			throw "File: $arg does not exist\n";
		}

		push( @{ $self->{known} }, $arg );
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


1;
