package ReseqTrack::DBSQL::ERADBAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::DBConnection;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::DBSQL::ERARunMetaInfoAdaptor;

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

1;
