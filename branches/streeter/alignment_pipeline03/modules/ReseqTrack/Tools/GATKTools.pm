
=pod

=head1 NAME

ReseqTrack::Tools::GATKTools

=head1 SYNOPSIS

This is a base class for GATK tools e.g. IndelRealigner and QualityScoreRecalibrator

Takes a bam file as input and makes a bam file as output

ISA = ReseqTrack::Tools::RunProgram

If ENV's $SAMTOOLS and $GATK are set, then no need to 
pass to object. They will be assigned accordingly.

SAMTOOLS required to index bams and extract headers.

example

my $GATK = $IR->new(
     -java_exe        =>"/usr/bin/java" ,
     -jvm_args        =>"-Xmx4g",
     -GATK_PATH       =>$GATK/GenomeAnalysisTK/",
     -working_dir     => '/path/to/dir/',
     -reference       => '/path/to/reference',
     -known_sites_files => '/path/to/file',
);



=cut


package ReseqTrack::Tools::GATKTools;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::RunSamtools;
use File::Basename;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use vars qw(@ISA);

@ISA = qw(ReseqTrack::Tools::RunProgram);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);

	my ( $java_exe, $jvm_args, $options, $gatk_path ,
	     $reference, $samtools, $known_sites_files) = rearrange(
		[
			qw(
			  JAVA_EXE
			  JVM_ARGS
			  OPTIONS
			  GATK_PATH
                          REFERENCE
                          SAMTOOLS
                          KNOWN_SITES_FILES
			  )
		],
		@args
	);

	#defaults
        if (!$self->program && !$java_exe) {
          $java_exe = 'java';
        }
	$jvm_args ||= "-Xmx4g";
        $gatk_path ||= $ENV{GATK};

	$self->java_exe($java_exe);
	$self->jvm_args($jvm_args);
	$self->gatk_path($gatk_path);
	$self->reference($reference);
        $self->samtools($samtools);
        $self->known_sites_files($known_sites_files);

        while (my @key_value = each(%$options)) {
          $self->options(@key_value);
        }

	return $self;
}


sub generate_job_name {
    my $self = shift;

    my $job_name = basename($self->input_bam);
    $job_name =~ s/(\w+).*/$1/;
    $job_name .= $$;

    $self->job_name($job_name);
    return $job_name;
}

sub check_bai_exists{
  my ($self) = @_;

  my $bamindex = $self->input_bam . "\.bai";
  return if (-e $bamindex);

  print "$bamindex does not exist. Creating\n";

  my $samtools_object = ReseqTrack::Tools::RunSamtools->new(
                -program => $self->samtools,
                -input_files => $self->input_bam,
                        );
  $samtools_object->run('index');
  $bamindex = $samtools_object->output_files->[0];
  $self->created_files($bamindex);

  print "Created $bamindex\n\n";
  return;

}





sub input_bam {
  my $self = shift;
  return $self->input_files(@_)->[0];
}

sub output_bam {
  my $self = shift;
  return $self->output_files(@_)->[0];
}

sub check_jar_file_exists {
	my ($self) = shift;
	my $jar_file = $self->gatk_path . "\/" . $self->jar_file;
        check_file_exists($jar_file);
	return;
}

sub java_exe {
  my $self = shift;
  return $self->program(@_);
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

sub samtools {
  my ( $self, $arg ) = @_;
  if ($arg) {
    check_file_exists($arg);
    $self->{samtools} = $arg;
  }
  return $self->{samtools};
}

sub reference {
  my ( $self, $arg ) = @_;
  if ($arg) {
    check_file_exists($arg);
    $self->{reference} = $arg;
  }
  return $self->{reference};

}

sub known_sites_files {
  my ( $self, $arg ) = @_;

  $self->{'known_sites_files'} ||= {};
  if ($arg) {
    foreach my $file (@{ref($arg) eq 'ARRAY' ? $arg : [$arg]}) {
      $file =~ s{//}{/}g;
      $self->{'known_sites_files'}->{$file} = 1;
    }
  }

  my @files = keys %{$self->{'known_sites_files'}};
  return \@files;
}

1;
