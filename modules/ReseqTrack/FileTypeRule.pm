
=pod

=head1 NAME

ReseqTrack::FileTypeRule

=head1 SYNOPSIS

This is a container object for the file_type_rule table.
The file_type_rule table describes the rules that are applied to assign a file type, and the order in which they are applied

=head1 Example

my $file_type_rule = ReseqTrack::FileTypeRule->new(
      -rule_block_order => $rule_block_order,
      -rule_order => $rule_order,
      -file_type  => $file_type,
      -match_regex  => $match_regex,
        );

=cut

package ReseqTrack::FileTypeRule;

use strict;
use warnings;

use base qw(ReseqTrack::Base);

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);


=head2 new
  
  Arg [1]    : ReseqTrack::FileTypeRule
  Arg [2]    : int, order in which a block of rules must be considered
  Arg [3]    : int, order in which the rule must be considered with the block
  Arg [4]    : string, file type
  Arg [5]    : string, regular expressions for the filename to match, separated by brackets and keywords NOT, AND, OR
  Function   : create a ReseqTrack::FileTypeRule object
  Returntype : ReseqTrack::FileTypeRule
  Exceptions : throws if no rule_block_order or rule_order specified
  Example    : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my (
        $rule_block_order,  $rule_order,     $file_type,      $match_regex,
      )
      = rearrange(
        [
            qw(RULE_BLOCK_ORDER RULE_ORDER FILE_TYPE MATCH_REGEX)
        ],
        @args
      );

    #ERROR CHECKING
    throw("Can't create ReseqTrack::FileTypeRule without a rule_block_order") unless ($rule_block_order);
    throw("Can't create ReseqTrack::FileTypeRule without a rule_order") unless ($rule_order);
    ######
    $self->rule_block_order($rule_block_order);
    $self->rule_order($rule_order);
    $self->file_type($file_type);
    $self->match_regex($match_regex);
    #########
    return $self;
}



=head2 accessor methods

  Arg [1]   : ReseqTrack::FileTypeRule
  Arg [2]   : variable
  Function  : store variable in object
  Returntype: variable
  Exceptions: none
  Example   : 

=cut


sub rule_block_order {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'rule_block_order'} = $arg;
    }
    return $self->{'rule_block_order'};
}

sub rule_order {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'rule_order'} = $arg;
    }
    return $self->{'rule_order'};
}

sub file_type {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'file_type'} = $arg;
    }
    return $self->{'file_type'};
}

sub match_regex {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{'match_regex'} = $arg;
        if ($self->{'perl_condition'}) {
            delete $self->{'perl_condition'};
        }
    }
    return $self->{'match_regex'};
}

=head2 perl_condition

  Arg [1]   : ReseqTrack::FileTypeRule
  Arg [2]   : string, perl variable name, defaults to 'filename'
  Function  : gets the filename match condition in the form of perl code. Checks if match_regex has already been converted to perl code.
  Returntype: string, perl code
  Exceptions: none
  Example   : my $perl_condition = $file_type_rule->perl_condition('filename');

=cut
sub perl_condition {
    my ( $self, $perl_variable_name ) = @_;
    $perl_variable_name = $perl_variable_name || 'filename';

    if (! $self->{'perl_condition'}{$perl_variable_name}) {
        $self->{'perl_condition'}{$perl_variable_name}
                = $self->create_perl_condition($perl_variable_name);
    }
    return $self->{'perl_condition'}{$perl_variable_name};
}

sub create_perl_condition {
    my ($self, $perl_variable_name) = @_;
    return 1 if (!$self->match_regex);

    my $match_regex = $self->match_regex;
    $match_regex =~ s/NOT\s+(\/\S*\/i?)/NOT \($1\)/ig;

    my @split_match_regex = split(/(\/\S*\/i?)/, $match_regex);
    foreach my $match_part (@split_match_regex) {
        if ($match_part =~ /\/\S*\/i?/) {
            $match_part = '$' . $perl_variable_name . ' =~ ' . $match_part;
        }
        else {
            $match_part =~ s/AND/\&\&/ig;
            $match_part =~ s/OR/\|\|/ig;
            $match_part =~ s/NOT/!/g;
        }
    }
    my $perl_condition = join('', @split_match_regex);
    return $perl_condition;
}

=head2 test_filename

  Arg [1]   : ReseqTrack::FileTypeRule
  Arg [2]   : string, filename
  Function  : Tests whether the match_regex matches the filename
  Returntype: 1 / 0
  Example   : my $matched = $file_type_rule->test_filename($filename);

=cut

sub test_filename {
    my ($self, $filename) = @_;
    my $perl_condition = $self->perl_condition('filename');

    my $matched = eval "$perl_condition";
    return $matched;
}







1;
