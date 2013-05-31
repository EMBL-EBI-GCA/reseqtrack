=pod

=head1 NAME

ReseqTrack::HasAttributes

=head1 SYNOPSIS

Base Class for objects which can have attribute objects attached to them

=head1 Example

my $file = ReseqTrack::Sample->new(
      -name => $path,
      -type => $type,
      -size => $size,
      -host => $host,
      -attributes => \%attributes,
        );


=cut
package ReseqTrack::HasAttributes;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);

use ReseqTrack::Base;

@ISA = qw(ReseqTrack::Base);

=head2 new

  Arg [1]   : ReseqTrack::HasAttributes
  Arg [2]   : hashref of key value pairs
  Function  : create ReseqTrack::HasAttributes object
  Returntype: ReseqTrack::HasAttributes
  Exceptions: 
  Example   : 

=cut

sub new {
    my ( $class, @args ) = @_;
    my $self = $class->SUPER::new(@args);
    my ( $attributes ) =
      rearrange( [ 'ATTRIBUTES'], @args );
    $self->attributes($attributes) if ($attributes);
    return $self;
}

sub attribute{
	my ($self,$name, $value) = @_;
	
	if (!defined $name ) {
		throw("Must pass an attribute name");
	}
	if (defined $value){
		$self->{attributes}->{$name} = $value;
	}
	
	return $self->{attributes}->{$name};
}

sub attributes{
	my ($self,$attributes) = @_;
	if (ref($attributes) eq 'HASH') {
		$self->{attributes} = $attributes;
	}
	else {
		throw("Arg to attributes() must be a hash ref")
	}
}