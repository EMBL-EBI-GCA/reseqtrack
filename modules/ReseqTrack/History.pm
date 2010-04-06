
=pod

=head1 NAME

ReseqTrack::History;

=head1 SYNOPSIS

An container object to hold the data for the history table.

The history table is designed to store information about changes to other
objects in other tables like movement of files or updates to run meta info data

my $history = ReseqTrack::History->new(
                -table_name => 'file',
                -comment => 'moved from /one/path to /two/path',
                 );

=cut

package ReseqTrack::History;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

=head2 new

  Arg [1]   : int, id
  Arg [2]   : string, table name
  Arg [3]   : string, comment
  Arg [4]   : string, date time stamp
  Function  : create ReseqTrack::History object
  Returntype: ReseqTrack::History
  Exceptions: throws if not given an other id, a table name or if the comment
  is an empty string
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my ( $other_id, $table_name, $comment, $time ) =
      rearrange( [qw(OTHER_ID TABLE_NAME COMMENT TIME)], @args );

    #error checking
    throw("Can't create history object wihtout another id") unless ($other_id);
    throw("Can't create a history object with an undefined table name")
      unless ($table_name);
    throw( "Can't create a history object with a table name " . $table_name )
      unless ( $table_name eq 'file'
        || $table_name eq 'event'
        || $table_name eq 'run_meta_info'
        || $table_name eq 'alignment_meta_info'
        || $table_name eq 'collection' );
    throw("Can't create a history object without a comment")
      unless ($comment);
    ######
    $self->other_id($other_id);
    $self->table_name($table_name);
    $self->comment($comment);
    $self->time($time);

    return $self;
}

=head2 accessor methods

  Arg [1]   : Reseq::History
  Arg [2]   : mostly strings
  Function  : set variable in object
  Returntype: return variable
  Exceptions: n/a
  Example   : 

=cut

sub other_id {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{other_id} = $arg;
    }
    return $self->{other_id};
}

sub table_name {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{table_name} = $arg;
    }
    return $self->{table_name};
}

sub comment {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{comment} = $arg;
    }
    return $self->{comment};
}

sub time {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{time} = $arg;
    }
    return $self->{time};
}

sub object_table_name {
    my ($self) = @_;
    return 'history';
}

1;
