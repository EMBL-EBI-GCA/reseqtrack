package ReseqTrack::DBSQL::ERADBAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::DBConnection;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::DBSQL::ERARunMetaInfoAdaptor;
use ReseqTrack::DBSQL::ERA::StudyAdaptor;
use ReseqTrack::DBSQL::ERA::RunAdaptor;
use ReseqTrack::DBSQL::ERA::ExperimentAdaptor;
use ReseqTrack::DBSQL::ERA::SampleAdaptor;
use base qw(ReseqTrack::DBSQL::DBAdaptor);
use Data::Dumper;

sub new{
  # print "test_again"
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
  print Dumper($self);
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
	return $self->_get_adaptor('ReseqTrack::DBSQL::ERA::StudyAdaptor');
}
sub get_SampleAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::ERA::SampleAdaptor');
}
sub get_ExperimentAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::ERA::ExperimentAdaptor');
}
sub get_RunAdaptor{
	my ($self) = @_;
	return $self->_get_adaptor('ReseqTrack::DBSQL::ERA::RunAdaptor');
}


1;
