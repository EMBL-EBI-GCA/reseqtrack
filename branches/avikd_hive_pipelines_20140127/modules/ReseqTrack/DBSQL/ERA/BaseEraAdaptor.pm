package ReseqTrack::DBSQL::ERA::BaseEraAdaptor;

use strict;
use warnings;
use base qw(ReseqTrack::DBSQL::BaseAdaptor);
use ReseqTrack::Tools::Exception qw(throw);
use XML::Twig;

use ReseqTrack::Tools::AttributeUtils qw(create_attribute_for_object);

sub store {
  throw("We do not update ERAPRO");
}

sub update {
  throw("We do not update ERAPRO");
}

sub fetch_by_source_id {
  my ( $self, $source_id ) = @_;
  return $self->fetch_by_dbID($source_id);
}

sub xml_column {
  throw("Specify xml column name");
}

sub attach_attributes {
  my ( $self, $object ) = @_;
  my $xml_column         = $self->xml_column();
  my $internal_id_column = $self->internal_id_column();
  my $source_id          = $object->source_id;

  my $sql =
    "select  xmltype.getclobval($xml_column) xml from " . $self->table_name;
  $sql .= " where " . $self->where;
  $sql .= " and $internal_id_column = ?";

  my $sth = $self->prepare($sql);
  $sth->bind_param( 1, $source_id );
  eval { $sth->execute; };
  if ($@) {
    throw("Problem running $sql $@");
  }

  my $ary_ref = $sth->fetch;
  throw( "No row returned with $sql and " . $source_id )
    unless $ary_ref;

  $self->attach_attributes_from_xml( $object, $ary_ref->[0] );

  $sth->finish;
}

sub attach_attributes_from_xml {
  my ( $self, $object, $xml ) = @_;

  my $attribute_tag = $self->attribute_tag;
  my @attributes;

  my $twig = XML::Twig->new(
    twig_handlers => {
      $attribute_tag => sub {
        my ( $t, $element ) = @_;

        my $key   = uc( $element->first_child_text('TAG') );
        my $value = $element->first_child_text('VALUE');
        
        for ( $key, $value ) {
          s/^\s+|\s+$//g;
        }

        if ($value) {
          push @attributes,
            create_attribute_for_object( $object, $key, $value );
        }
       }
    }
  );

  $twig->parse($xml);
  $object->attributes( \@attributes );
}

sub fetch_by_study_id {
  throw("Implement me");
}

sub add_ega_id {
  my ( $self, $object, $hashref ) = @_;

  if ( $hashref->{EGA_ID} ) {
    my $attr =
      create_attribute_for_object( $object, 'EGA_ID', $hashref->{EGA_ID} );
    $object->attributes( [$attr] );
  }
}

1;
