package ReseqTrack::GenotypeResults;

use strict;
use warnings;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning);
use vars qw(@ISA);

@ISA = qw(ReseqTrack::Base);

sub new {

 my ( $class, @args ) = @_;
 my $self = $class->SUPER::new(@args);

  my (

      $genotype_results_id,
      $other_id,
      $table_name,
      $reference,
      $snps_bin,
      $aligner,
      $version,
      $validation_method,
      $max_bases,
      $percent_mapped,
      $summary,
      $verdict,
      $performed,

     ) = rearrange(
		   [
		    qw(
			  GENOTYPE_RESULTS_ID
			  OTHER_ID
			  TABLE_NAME
			  REFERENCE
			  SNPS_BIN
			  ALIGNER
			  VERSION
			  VALIDATION_METHOD
			  MAX_BASES
			  PERCENT_MAPPED
			  SUMMARY
			  VERDICT
			  PERFORMED )
		   ],
		   @args
		  );
 
  ######
  $self->genotype_results_id($genotype_results_id);
  $self->other_id($other_id);
  $self->table_name($table_name);
  $self->reference($reference);
  $self->snps_bin($snps_bin);
  $self->aligner($aligner);
  $self->version($version);
  $self->validation_method($validation_method);
  $self->max_bases($max_bases);
  $self->percent_mapped($percent_mapped);
  $self->summary($summary);
  $self->verdict($verdict);
  $self->performed($performed);

  #########

  my $ERR_MSG = "Can't create ReseqTrack::Genotyping without";

  #ERROR CHECKING
  throw("$ERR_MSG other_id")   unless ($other_id);
  throw("$ERR_MSG table_name") unless ($table_name);
  throw("$ERR_MSG reference")  unless ($reference);
  throw("$ERR_MSG snps_bin")   unless ($snps_bin);
  throw("$ERR_MSG aligner")    unless ($aligner);
  throw("$ERR_MSG version")    unless ($version);
  throw("$ERR_MSG validation_method")    unless ($validation_method);
  throw("$ERR_MSG summary")    unless ($summary);
  #throw("$ERR_MSG max_bases")  unless ($max_bases);
  #throw("$ERR_MSG verdict")    unless ($verdict);

  return $self;
}

=head2 accessor methods

  Arg [1]   : ReseqTrack::Genotyping
  Arg [2]   : variable
  Function  : store variable in object
  Returntype: variable
  Exceptions: none
  Example   : 

=cut
sub genotype_results_id {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{genotype_results_id} = $arg;
  }
  return $self->{genotye_results_id};
}
sub other_id {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{other_id} = $arg;
  }
  return $self->{other_id};
}

sub table_name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{table_name} = $arg;
  }
  return $self->{table_name};
}

sub aligner {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{aligner} = $arg;
  }
  return $self->{aligner};
}
sub version {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{version} = $arg;
  }
  return $self->{version};
}

sub validation_method {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{validation_method} = $arg;
  }
  return $self->{validation_method};
}

sub max_bases {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{max_bases} = $arg;
  }
  return $self->{max_bases};
}

sub percent_mapped {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{percent_mapped} = $arg;
  }
  return $self->{percent_mapped};
}
sub performed {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{preformed} = $arg;
  }
  return $self->{performed};
}

sub snps_bin {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{snps_bin} = $arg;
  }
  return $self->{snps_bin};
}

sub verdict {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{verdict} = $arg;
  }
  return $self->{verdict};
}

sub reference {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{reference} = $arg;
  }
  return $self->{reference};
}

sub summary {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{summary} = $arg;
  }
  return $self->{summary};
}

1;
