package ReseqTrack::DBSQL::FileTypeRuleAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::FileTypeRule;
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
  return "file_type_rule.rule_block_order, file_type_rule.rule_order, file_type_rule.file_type, file_type_rule.match_regex";
}

sub table_name{
  return "file_type_rule";
}


sub store{
  my ($self, $file_type_rule) = @_;
  my $sql = "insert ignore into file_type_rule ".
      "(rule_block_order, rule_order, file_type, match_regex) ".
      "values(?, ?, ?, ?) ";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $file_type_rule->rule_block_order);
  $sth->bind_param(2, $file_type_rule->rule_order);
  $sth->bind_param(3, $file_type_rule->file_type);
  $sth->bind_param(4, $file_type_rule->match_regex);
  my $rows_inserted = $sth->execute();
  $sth->finish();
  $file_type_rule->adaptor($self);

  delete $ARRAY_CACHE{$self->db};
}



sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create a ReseqTrack::FileTypeRule from an undefined hashref") if(!$hashref);
  my $file_type_rule = ReseqTrack::FileTypeRule->new(
				   -rule_block_order => $hashref->{rule_block_order},
				   -rule_order => $hashref->{rule_order},
				   -file_type => $hashref->{file_type},
				   -match_regex => $hashref->{match_regex},
				   -adaptor => $self,
				  );
  return $file_type_rule;
}

sub fetch_all_in_order{
    my ($self, $no_cache) = @_;

    if (!$no_cache && $ARRAY_CACHE{$self->db}) {
        return $ARRAY_CACHE{$self->db};
    }

    my $file_type_rules = $self->fetch_all;

    my %ftr_hash;
    foreach my $ftr (@$file_type_rules) {
        $ftr_hash{$ftr->rule_block_order}{$ftr->rule_order} = $ftr;
    }

    my @ordered_file_type_rules;
    foreach my $block_order (sort {$a <=> $b} keys %ftr_hash) {
        my @ordered_block = map {$ftr_hash{$block_order}{$_}} sort {$a <=> $b} keys %{$ftr_hash{$block_order}};
        push(@ordered_file_type_rules, \@ordered_block);
    }

    if (!$no_cache) {
        $ARRAY_CACHE{$self->db} = \@ordered_file_type_rules;
    }
    return \@ordered_file_type_rules;
}

sub delete_all{
  my $self = shift;
  my $sql = "delete from file_type_rule";
  my $sth = $self->prepare($sql);
  $sth->execute();
  $sth->finish();
  delete $ARRAY_CACHE{$self->db};
}


1;
