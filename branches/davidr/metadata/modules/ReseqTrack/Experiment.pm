package ReseqTrack::Experiment;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::HasHistory;

@ISA = qw(ReseqTrack::HasHistory);


sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  my ($study_id, $status, $md5, $center_name, 
  	$experiment_alias, $instr_platform, $instr_model, $lib_layout, 
  	$lib_name, $lib_strategy, $lib_source, $lib_selection, 
  	$paired_nominal_length, $paired_nominal_sdev, $source_id ) =
      rearrange([qw( STUDY_ID STATUS MD5 CENTER_NAME
    EXPERIMENT_ALIAS INSTRUMENT_PLATFORM INSTRUMENT_MODEL LIBRARY_LAYOUT     
    LIBRARY_NAME LIBRARY_STRATEGY LIBRARY_SOURCE LIBRARY_SELECTION  
    PAIRED_NOMINAL_LENGTH PAIRED_NOMINAL_SDEV SOURCE_ID) ], @args);
  
	$self->study_id($study_id);
 	$self->status($status);
	$self->md5($md5);
	$self->center_name($center_name); 
	$self->experiment_alias($experiment_alias);
  $self->library_name($lib_name);
  $self->library_strategy($lib_strategy);
  $self->library_source($lib_source);
  $self->library_selection($lib_selection); 
  $self->paired_nominal_length($paired_nominal_length);
  $self->instrument_platform($instr_platform);
	$self->instrument_model($instr_model);
	$self->library_layout($lib_layout); 
  $self->paired_nominal_sdev($paired_nominal_sdev);
  $self->source_id($source_id); 
  
  return $self;
}

sub source_id{
  my ($self, $arg) = @_; 
  
  if($arg){
    $self->{source_id} = $arg;
  }
  return $self->{source_id};
}

sub study{
	my ($self, $arg) = @_;
	
	if ($arg){
		$self->{study} = $arg;
	}
	elsif ($self->adaptor() && $self->study_id()){
		my $sa = $self->adaptor()->get_StudyAdaptor();
		$self->{study} = $sa->fetch_by_dbID( $self->study_id() ); 
	}
	
	return $self->{study};
}

sub study_id{
  my ($self, $arg) = @_;
  if($arg){
    $self->{study_id} = $arg;
  }
  return $self->{study_id};
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

sub experiment_alias{
  my ($self, $arg) = @_;
  if($arg){
    $self->{experiment_alias} = $arg;
  }
  return $self->{experiment_alias};
}

sub instrument_platform{
  my ($self, $arg) = @_;
  if($arg){
    $self->{instrument_platform} = $arg;
  }
  return $self->{instrument_platform};
}

sub instrument_model{
  my ($self, $arg) = @_;
  if($arg){
    $self->{instrument_model} = $arg;
  }
  return $self->{instrument_model};
}   

sub library_layout{
  my ($self, $arg) = @_;
  if($arg){
    $self->{library_layout} = $arg;
  }
  return $self->{library_layout};
}     

sub library_name{
  my ($self, $arg) = @_;
  if($arg){
    $self->{library_name} = $arg;
  }
  return $self->{library_name};
}      
 
sub library_strategy{
  my ($self, $arg) = @_;
  if($arg){
    $self->{library_strategy} = $arg;
  }
  return $self->{library_strategy};
}  
 
sub library_source{
  my ($self, $arg) = @_;
  if($arg){
    $self->{library_source} = $arg;
  }
  return $self->{library_source};
}  
   
sub library_selection{
  my ($self, $arg) = @_;
  if($arg){
    $self->{library_selection} = $arg;
  }
  return $self->{library_selection};
}  

sub paired_nominal_length{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{paired_nominal_length} = $arg;
  }
  return $self->{paired_nominal_length};
}

sub paired_nominal_sdev{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{paired_nominal_sdev} = $arg;
  }
  return $self->{paired_nominal_sdev};
}

sub object_table_name {
	return "experiment";
}

1;
