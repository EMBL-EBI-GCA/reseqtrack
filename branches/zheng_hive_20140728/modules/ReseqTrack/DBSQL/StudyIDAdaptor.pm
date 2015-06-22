package ReseqTrack::DBSQL::StudyIDAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::Tools::Exception qw(throw warning);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub columns{
  return "study_id.study_id";
}

sub table_name{
  return "study_id";
}


sub store{
  my ($self, $study_id) = @_;
  my $sql = "insert ignore into study_id ".
      "(study_id) ".
      "values(?) ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $study_id);
  my $rows_inserted = $sth->execute();
  $sth->finish();
}


sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create a study_id from an undefined hashref") if(!$hashref);
  return $hashref->{study_id};
}

1;
