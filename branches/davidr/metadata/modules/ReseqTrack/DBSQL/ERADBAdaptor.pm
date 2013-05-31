package ReseqTrack::DBSQL::ERADBAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::DBConnection;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::DBSQL::ERARunMetaInfoAdaptor;
use ReseqTrack::DBSQL::StudyAdaptor;
use ReseqTrack::DBSQL::RunAdaptor;
use ReseqTrack::DBSQL::ExperimentAdaptor;
use ReseqTrack::DBSQL::SampleAdaptor;

@ISA = qw(ReseqTrack::DBSQL::DBAdaptor);

sub new{
  my ($class, @args) = @_;
  my $has_driver = 0;
  foreach my $arg(@args){
    if($arg && $arg =~ /driver/i){
      $has_driver = 1;
    }
  }
  unless($has_driver){
    push(@args, ('-driver', 'Oracle'));
  }
  my $self = $class->SUPER::new(@args);
  return $self;
}

sub get_ERARunMetaInfoAdaptor{
  my ($self, @adaptor_args) = @_;
  my ($study_ids, $population_rules) =
    rearrange([qw(STUDY_IDS POPULATION_RULES)], @adaptor_args);
  if($self->{run_meta_info_adaptor}){
      if ($study_ids) {
          $self->{run_meta_info_adaptor}->study_ids($study_ids);
      }
      if ($population_rules) {
          $self->{run_meta_info_adaptor}->study_ids($study_ids);
      }
  }
  else {
    $self->{run_meta_info_adaptor} = ReseqTrack::DBSQL::ERARunMetaInfoAdaptor->
        new(-db => $self, @adaptor_args);
  }
  return $self->{run_meta_info_adaptor};
}

sub get_StudyAdaptor{
  my ($self) = @_;
  if(!$self->{study_adaptor}){
    $self->{study_adaptor} = ReseqTrack::DBSQL::StudyAdaptor->new($self);
  }
  return $self->{study_adaptor};
}

sub get_SampleAdaptor{
  my ($self) = @_;
  if(!$self->{sample_adaptor}){
    $self->{sample_adaptor} = ReseqTrack::DBSQL::SampleAdaptor->new($self);
  }
  return $self->{sample_adaptor};
}

sub get_RunAdaptor{
	my ($self) = @_;
  if(!$self->{run_adaptor}){
    $self->{run_adaptor} = ReseqTrack::DBSQL::RunAdaptor->new($self);
  }
  return $self->{run_adaptor};
}

sub get_ExperimentAdaptor{
	my ($self) = @_;
  if(!$self->{experiment_adaptor}){
    $self->{experiment_adaptor} = ReseqTrack::DBSQL::ExperimentAdaptor->new($self);
  }
  return $self->{experiment_adaptor};
}

1;
