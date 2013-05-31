package ReseqTrack::Sample;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use Data::Dumper;
use ReseqTrack::HasAttributes;

@ISA = qw(ReseqTrack::HasAttributes);

sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($status, $md5, $center_name, 
		$sample_alias, $tax_id, $scientific_name, $common_name, 
		$anonymized_name, $individual_name, $sample_title,) =
      rearrange([qw(STATUS MD5 CENTER_NAME
      SAMPLE_ALIAS TAX_ID SCIENTIFIC_NAME COMMON_NAME
			ANONYMIZED_NAME INDIVIDUAL_NAME SAMPLE_TITLE) ], @args);

  $self->status($status);
	$self->md5($md5);
	$self->center_name($center_name); 
	$self->sample_alias($sample_alias);
	$self->tax_id($tax_id);
	$self->scientific_name($scientific_name);
	$self->common_name($common_name); 
	$self->anonymized_name($anonymized_name);
	$self->individual_name($individual_name);
	$self->sample_title($sample_title);
  
  return $self;
}
 
sub status{
  my ($self, $arg) = @_;
  if($arg){
    $self->{status} = $arg;
  }
  return $self->{status};
}

sub md5{
  my ($self, $arg) = @_;
  if($arg){
    $self->{md5} = $arg;
  }
  return $self->{md5};
}

sub center_name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{center_name} = $arg;
  }
  return $self->{center_name};
}

sub sample_alias{
  my ($self, $arg) = @_;
  if($arg){
    $self->{sample_alias} = $arg;
  }
  return $self->{sample_alias};
}

sub tax_id{
  my ($self, $arg) = @_;
  if($arg){
    $self->{tax_id} = $arg;
  }
  return $self->{tax_id};
}

sub scientific_name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{scientific_name} = $arg;
  }
  return $self->{scientific_name};
}

sub common_name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{common_name} = $arg;
  }
  return $self->{common_name};
}

sub anonymized_name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{anonymized_name} = $arg;
  }
  return $self->{anonymized_name};
}
 
sub individual_name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{individual_name} = $arg;
  }
  return $self->{individual_name};
}

sub sample_title{
  my ($self, $arg) = @_;
  if($arg){
    $self->{sample_title} = $arg;
  }
  return $self->{sample_title};
}

1;
