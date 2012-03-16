=pod

=head1 NAME

ReseqTrack::Tools::GATKTools::IndelReAligner

=head1 SYNOPSIS

Object to create a bam file that in realigned around
known indel sites using GATK RealignerTargetCreator and
IndelRealigner

example

my $REALIGN_AROUND_INDELS = $IR->new(
		     -reference       => $reference,
		     -input_files     => $input{input_files},
		     -rtc_knowns      =>" -known $millsindels.vcf.gz ",
		     -working_dir     => $input{working_dir},
);



=cut

package ReseqTrack::Tools::GATKTools::IndelReAligner;

use strict;
use warnings;

use ReseqTrack::Tools::RunProgram;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use ReseqTrack::Tools::GATKTools;

use vars qw(@ISA);
use Data::Dumper;

@ISA = qw(ReseqTrack::Tools::GATKTools);

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);


        #setting defaults
        if (!$self->jar_file) {
          $self->jar_file ("GenomeAnalysisTK.jar");
        }

        if (!$self->options('IndelRealigner')) {
          my $ira_options = " -LOD 0.4 -model KNOWNS_ONLY --disable_bam_indexing";
          $self->options('IndelRealigner', $ira_options);
        }

	return $self;
}


sub run_program {
	my $self = shift;

        throw "no input bam" if (!$self->input_bam);
        throw "\nNo IndelRealigner options\n"
            if ( ! $self->{options}->{'IndelRealigner'});
        warn "No known indels files"
            if (! @{$self->known_sites_files});
        check_file_exists($_) foreach (@{$self->known_sites_files});
        $self->check_jar_file_exists;
        check_file_exists($self->reference);
        $self->check_bai_exists();


	$self->create_target_intervals_file();
	$self->create_indel_realign_bam();

	return;
}

sub check_bai_exists{
  my ($self) = @_;

  my $bamindex = $self->input_bam . "\.bai";
  return if (-e $bamindex);

  print "$bamindex does not exist. Creating\n";

  my $samtools_object = ReseqTrack::Tools::RunSamtools->new(
                -program => $self->samtools, -flag_index => 1,
                -input_files => $self->input_bam,
                        );
  $samtools_object->run;
  $bamindex = $samtools_object->output_bai_files->[0];
  $self->created_files($bamindex);

  print "Created $bamindex\n\n";
  return;

}


sub create_target_intervals_file {
  my $self = shift;

  my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
  $cmd .= $self->gatk_path ;
  $cmd .= "\/";
  $cmd .= $self->jar_file;
  $cmd .= " -T RealignerTargetCreator ";

  if ( defined 	$self->{options}->{'RealignerTargetCreator'}) {
    $cmd .= 	$self->{options}->{'RealignerTargetCreator'};
  }
  else {
    warn "No RealignerTargetCreator options\n";
  }

  foreach my $vcf (@{$self->known_sites_files}) {
    $cmd .= "-known $vcf ";
  }

  my $interval_file = $self->working_dir . '/'
        . $self->job_name . '.bam.interval_list';
  $interval_file =~ s{//}{/}g;

  $cmd .= " -o $interval_file ";
  $cmd .= " -R " . $self->reference . " ";
  $cmd .= " -I " . $self->input_bam;

  my $CL = "java jvm_args -jar GenomeAnalysisTK.jar "
    . "-T RealignerTargetCreator -R \$reference -o \$intervals_file "
      ."-known \$known_sites_file(s)";

  #For bam header section
  my %PG =('PG'=>'@PG',
           'ID'=>"gatk_target_interval_creator",
           'PN'=>"GenomeAnalysisTK",     
           'PP'=>"sam_to_fixed_bam",
           'VN'=>"1.2-29-g0acaf2d",
           'CL'=> $CL,
           );

  my $COMMENT = '$known_sites_file(s) = '. join(', ', @{$self->known_sites_files});
  my %CO = ('@CO'=> $COMMENT);

  $self->intervals_file($interval_file);
  $self->execute_command_line ($cmd);


  return;
}


sub create_indel_realign_bam {
  my $self = shift;

  my $cmd = $self->java_exe . " " . $self->jvm_args . " -jar ";
  $cmd .= $self->gatk_path . "\/" . $self->jar_file;
  $cmd .= " -T IndelRealigner ";
  $cmd .= $self->{options}->{'IndelRealigner'};
  $cmd .= " --targetIntervals " . $self->intervals_file . " ";
  foreach my $vcf (@{$self->known_sites_files}) {
    $cmd .= "-known $vcf ";
  }
  $cmd .= " -R " . $self->reference . " ";
  $cmd .= " -I " . $self->input_bam;

  my $realigned_bam = $self->working_dir . '/'
        . $self->job_name. '.indel_realigned.bam';
  $cmd .= " -o " . $realigned_bam;

  $self->output_files($realigned_bam);
  $self->execute_command_line ($cmd);

  return;

}

sub intervals_file {

  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{intervals_file} = $arg;
    $self->created_files($arg);
  }
  return $self->{intervals_file};
}

1;
