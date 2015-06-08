package ReseqTrack::DBSQL::PopulationRuleAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::PopulationRule;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

my %ARRAY_CACHE;

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub columns{
  return "population_rule.rule_order, population_rule.population, population_rule.match_regex";
}

sub table_name{
  return "population_rule";
}


sub store{
  my ($self, $population_rule) = @_;
  my $sql = "insert ignore into population_rule ".
      "(rule_order, population, match_regex) ".
      "values(?, ?, ?) ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $population_rule->rule_order);
  $sth->bind_param(2, $population_rule->population);
  $sth->bind_param(3, $population_rule->match_regex);
  my $rows_inserted = $sth->execute();
  $sth->finish();
  $population_rule->adaptor($self);

  delete $ARRAY_CACHE{$self->db};
}



sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create a ReseqTrack::PopulationRule from an undefined hashref") if(!$hashref);
  my $population_rule = ReseqTrack::PopulationRule->new(
				   -rule_order => $hashref->{rule_order},
				   -population => $hashref->{population},
				   -match_regex => $hashref->{match_regex},
				   -adaptor => $self,
				  );
  return $population_rule;
}

sub fetch_all_in_order{
    my ($self, $no_cache) = @_;

    if (!$no_cache && $ARRAY_CACHE{$self->db}) {
        return $ARRAY_CACHE{$self->db};
    }

    my $population_rules = $self->fetch_all;
    my @ordered_rules = sort {$a->rule_order <=> $b->rule_order} @$population_rules;

    if (!$no_cache) {
        $ARRAY_CACHE{$self->db} = \@ordered_rules;
    }
    return \@ordered_rules;
}

sub fetch_hash_of_populations{
  my ($self, $no_cache) = @_;
  my $rules = $self->fetch_all_in_order($no_cache);
  my %hash;
  foreach my $rule(@$rules){
    $hash{$rule->population} = 1;
  }
  return \%hash;
}

sub delete_all{
  my $self = shift;
  my $sql = "delete from population_rule";
  my $sth = $self->prepare($sql);
  $sth->execute();
  $sth->finish();
  delete $ARRAY_CACHE{$self->db};
}


1;
