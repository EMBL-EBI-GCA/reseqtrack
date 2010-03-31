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
  my ($self) = @_;
  if(!$self->{run_meta_info_adaptor}){
    $self->{run_meta_info_adaptor} = ReseqTrack::DBSQL::ERARunMetaInfoAdaptor->
        new($self);
  }
  return $self->{run_meta_info_adaptor};
}

1;
