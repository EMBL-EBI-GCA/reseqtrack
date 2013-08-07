
=pod

=head1 NAME

ReseqTrack::PopulationRule

=head1 SYNOPSIS

This is a container object for the population_rule table.
The population_rule table describes the rules that are applied to assign a population, and the order in which they are applied

=head1 Example

my $population_rule = ReseqTrack::PopulationRule->new(
      -rule_order => $rule_order,
      -population  => $population,
      -match_regex  => $match_regex,
        );

=cut

package ReseqTrack::PopulationRule;

use strict;
use warnings;

use base qw(ReseqTrack::Base);

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);


=head2 new

  Arg [1]    : ReseqTrack::PopulationRule
  Arg [2]    : int, order in which a rule must be considered
  Arg [3]    : string, population name
  Arg [4]    : string, regular expressions for a descriptive string to match, separated by brackets and keywords NOT, AND, OR
  Function   : create a ReseqTrack::PopulationRule object
  Returntype : ReseqTrack::PopulationRule
  Exceptions : throws if rule_order specified
  Example    : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my (
        $rule_order,     $population,      $match_regex,
      )
      = rearrange(
        [
            qw(RULE_ORDER POPULATION MATCH_REGEX)
        ],
        @args
      );

    #ERROR CHECKING
    throw("Can't create ReseqTrack::PopulationRule without a rule_order") unless ($rule_order);
    ######
    $self->rule_order($rule_order);
    $self->population($population);
    $self->match_regex($match_regex);
    #########
    return $self;
}



=head2 accessor methods

  Arg [1]   : ReseqTrack::PopulationRule
  Arg [2]   : variable
  Function  : store variable in object
  Returntype: variable
  Exceptions: none
  Example   : 

=cut

sub rule_order {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'rule_order'} = $arg;
    }
    return $self->{'rule_order'};
}

sub population {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'population'} = $arg;
    }
    return $self->{'population'};
}

sub match_regex {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'match_regex'} = $arg;
    }
    return $self->{'match_regex'};
}

=head2 perl_condition

  Arg [1]   : ReseqTrack::PopulationRule
  Function  : gets the description match condition in the form of perl code. Checks if match_regex has already been converted to perl code.
  Returntype: string, perl code
  Exceptions: none
  Example   : my $perl_condition = $population_rule->perl_condition();

=cut

sub perl_condition {
    my ( $self ) = @_;

    if (! $self->{'perl_condition'}) {
        if (!$self->match_regex) {
            $self->{'perl_condition'} = 1;
        }
        else {
            my $perl_condition = $self->match_regex;
            $perl_condition =~ s/\s+AND\s+/ \&\& /ig;
            $perl_condition =~ s/\s+OR\s+/ \|\| /ig;
            $perl_condition =~ s/\s+NOT\s+/ ! /g;
            $self->{'perl_condition'} = $perl_condition;
        }
    }
    return $self->{'perl_condition'};
}

=head2 test_description

  Arg [1]   : ReseqTrack::PopulationRule
  Arg [2]   : string, description
  Function  : Tests whether the match_regex matches the descriptive string
  Returntype: 1 / 0
  Example   : my $matched = $population_rule->test_description('han chinese');

=cut

sub test_description {
    my ($self, $description) = @_;
    my $perl_condition = $self->perl_condition();

    local $_ = $description; 
    my $matched = eval "$perl_condition";
    return $matched;
}







1;
