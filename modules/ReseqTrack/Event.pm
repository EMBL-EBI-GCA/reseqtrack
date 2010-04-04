
=pod

=head1 NAME

ReseqTrack::Event

=head1 SYNOPSIS

This is a container object for the event table. The event table describes program command lines
which can be run using the event pipeline system using input from several other tables in the
ReseqTrack database

=head1 Example

my $archive = ReseqTrack::Event->new
      (
       -name => 'archive_fastq_qc',
       -program  => 'perl run_fastq_qa.pl',
       -options => '-dbargs -output_dir /path/to/dir',
       -input_flag => '-run_id',
       -farm_options => '-q production',
       -output_path => '/path/to/output/',
       -type => 'FASTQ',  
       -table_name => 'file',
      );

=cut

package ReseqTrack::Event;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::HasHistory;

@ISA = qw(ReseqTrack::HasHistory);

=head2 new

  Arg [1]   : ReseqTrack::Event
  Arg [2]   : string, program name
  Arg [3]   : string, program version
  Arg [4]   : string, options for program
  Arg [5]   : string, command line flag for input string
  Arg [6]   : string, options for batch submission system
  Arg [7]   : int, batchsize this determines how many inputs will be run as part of the same
  farm job, this defaults to 0
  Arg [8]   : string, path to where output files should be written
  Arg [9]   : string, type the type string for the input from the given table
  Arg [10]   : string, table name to retrieve input from
  Arg [11]   : timestamp, time object was first stored in database
  Arg [12]   : timestamp, time object was last updated
  Function  : Create a ReseqTrack::Event object
  Returntype: ReseqTrack::Event
  Exceptions: throws if no name is defined
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);

    my (
        $name,       $program,      $program_version, $options,
        $input_flag, $farm_options, $batch_size,      $output_path,
        $type,       $table_name,   $created,         $updated
      )
      = rearrange(
        [
            qw(NAME
              PROGRAM
              PROGRAM_VERSION
              OPTIONS
              INPUT_FLAG
              FARM_OPTIONS
              BATCH_SIZE
              OUTPUT_PATH
              TYPE
              TABLE_NAME
              CREATED
              UPDATED)
        ],
        @args
      );

    throw("Can't create ReseqTrack::Event without a name") unless ($name);

    ######
    $self->name($name);                          #1
    $self->program($program);                    #2
    $self->program_version($program_version);    #3
    $self->options($options);                    #4
    $self->input_flag($input_flag);              #5
    $self->farm_options($farm_options);          #6
    $self->batch_size($batch_size);              #7
    $self->output_path($output_path);            #8
    $self->type($type);                          #9
    $self->table_name($table_name);              #10
    $self->created($created);                    #11
    $self->updated($updated);                    #12

    #########

    $self->batch_size(0) if ( !defined( $self->batch_size ) );

    return $self;
}

#It would be good to put in some docs here, look at ReseqTrack::File for the form

=head2 Accessor methods

  Arg [1]   : ReseaqTrack::Event
  Arg [2]   : string/int 
  Function  : set given arg to hash key of object
  Returntype: string/int
  Exceptions: n/a
  Example   : $event->name('archive_fastq_qc');

=cut

sub name {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{name} = $arg;
    }
    return $self->{name};
}

sub program {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{program} = $arg;
    }
    return $self->{program};
}

sub program_version {
    my ( $self, $arg ) = @_;

    if ($arg) {
        $self->{program_version} = $arg;
    }
    return $self->{program_version};
}

sub options {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{options} = $arg;
    }
    return $self->{options};
}

sub input_flag {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{input_flag} = $arg;
    }
    return $self->{input_flag};
}

sub farm_options {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{farm_options} = $arg;
    }
    return $self->{farm_options};
}

sub batch_size {
    my ( $self, $arg ) = @_;
    if ( defined($arg) ) {
        $self->{batch_size} = $arg;
    }
    return $self->{batch_size};
}

sub output_path {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{output_path} = $arg;
    }
    return $self->{output_path};
}

sub type {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{type} = $arg;
    }
    return $self->{type};
}

sub table_name {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{table_name} = $arg;
    }
    return $self->{table_name};
}

sub created {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{created} = $arg;
    }
    return $self->{created};
}

sub updated {
    my ( $self, $arg ) = @_;
    if ($arg) {
        $self->{updated} = $arg;
    }
    return $self->{updated};
}

=head2 object_table_name

  Arg [1]   : ReseqTrack::Event
  Function  : return table name for object, event
  Returntype: string
  Exceptions: 
  Example   : 

=cut

sub object_table_name {
    my ($self) = @_;
    return 'event';
}

1;
