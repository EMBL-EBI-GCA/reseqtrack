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

    
      $table_name, $other_id, $name,  $claimed,
      $top_hit, $second_hit, $ratio21,  $ratio_claimed,
      $reference, $snps_bin,  $aligner,  $version,
      $validation_method,  $max_bases, $percent_mapped,
      $percent_reads_used,  $verdict,  $cfg_file,
      $performed,

     ) = rearrange(
		   [
		    qw(
                  TABLE_NAME OTHER_ID NAME CLAIMED
                  TOP_HIT SECOND_HIT RATIO21 RATIO_CLAIMED
	          REFERENCE SNPS_BIN ALIGNER VERSION
		  VALIDATION_METHOD MAX_BASES PERCENT_MAPPED
                  PERCENT_READS_USED VERDICT CFG_FILE PERFORMED )
		   ],
		   @args
		  );
 

  $self->table_name($table_name);
  $self->other_id($other_id);
  $self->name($name);
  $self->claimed($claimed);
  $self->top_hit($top_hit);
  $self->second_hit($second_hit);
  $self->ratio21($ratio21);
  $self->ratio_claimed($ratio_claimed);
  $self->reference($reference);
  $self->snps_bin($snps_bin);
  $self->aligner($aligner);
  $self->version($version);
  $self->validation_method($validation_method);
  $self->max_bases($max_bases);
  $self->percent_mapped($percent_mapped);
  $self->percent_reads_used($percent_reads_used);
 
  $self->verdict($verdict);
  $self->performed($performed);
  $self->cfg_file($cfg_file);
  #########

  my $ERR_MSG = "Can't create ReseqTrack::GenotypeResults without";

  #ERROR CHECKING
  throw("$ERR_MSG other_id")   unless ($other_id);
  throw("$ERR_MSG table_name") unless ($table_name);
  throw("$ERR_MSG reference")  unless ($reference);
  throw("$ERR_MSG snps_bin")   unless ($snps_bin);
  throw("$ERR_MSG aligner")    unless ($aligner);
  throw("$ERR_MSG version")    unless ($version);
  throw("$ERR_MSG validation_method")    unless ($validation_method);
 
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
sub name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{name} = $arg;
  }
  return $self->{name};
}

sub claimed {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{claimed} = $arg;
  }
  return $self->{claimed};
}
sub top_hit {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{top_hit} = $arg;
  }
  return $self->{top_hit};
}
sub second_hit {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{second_hit} = $arg;
  }
  return $self->{second_hit};
}
sub ratio21 {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{ratio21} = $arg;
  }
  return $self->{ratio21};
}
sub ratio_claimed {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{ratio_claimed} = $arg;
  }
  return $self->{ratio_claimed};
}

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
sub percent_reads_used {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{percent_reads_used} = $arg;
  }
  return $self->{percent_reads_used};
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


sub cfg_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{cfg_file} = $arg;
  }
  return $self->{cfg_file};
}

1;