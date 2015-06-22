package ReseqTrack::Tools::Metadata::AddInBundle;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw);
use base qw(ReseqTrack::Tools::Metadata::BaseMetadataAddIn);

=head2 addins_to_load

  Function  : Provide a list of modules that should be loaded and wrapped by this bundle 
  Returntype: List of string, suggested new location for this file
	Implementation: See G1kBundle for an example implementation

=cut

sub addins_to_load {
  throw("Implement me!");
}

sub add_ins {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{add_ins} = $arg;
  }
  return $self->{add_ins};
}

sub new {
  my ( $class, @args ) = @_;

  my $this = $class->SUPER::new(@args);
  my $self = bless( $this, $class );

  my @add_ins;

  for my $module ( $self->addins_to_load ) {
    my $file = "$module.pm";
    $file =~ s{::}{/}g;
    eval { require "$file"; };
    if ($@) {
      throw("cannot load $file: $@");
    }

    my $add_in = $module->new(
      -era_db   => $self->era_db,
      -reseq_db => $self->reseq_db,
      -log_fh   => $self->log_fh
    );
    push @add_ins, $add_in;
  }

  $self->add_ins( \@add_ins );

  return $self;
}

sub check {
  my ( $self, $object, $last_copy ) = @_;

  for my $add_in ( @{ $self->add_ins } ) {
    $add_in->check( $object, $last_copy );
  }
}

sub report {
  my ($self) = @_;

  for my $add_in ( @{ $self->add_ins } ) {
    $add_in->report;
  }
}

sub isa {
  my ( $self, $class ) = @_;

  for my $add_in ( @{ $self->add_ins } ) {
    if ( $add_in->isa($class) ) {
      return 1;
    }
  }
  return 0;
}

1;
