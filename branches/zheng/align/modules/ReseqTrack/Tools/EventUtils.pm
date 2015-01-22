=pod

=head1 NAME

ReseqTrack::Tools::EventUtils;

=head1 SYNOPSIS

This is a collection of methods useful for running events

=head1 Example

use ReseqTrack::Tools::EventUtils qw (get_inputs);

my $input_hash = get_inputs($db, $events);

=cut

package ReseqTrack::Tools::EventUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::Job;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw();




1;
