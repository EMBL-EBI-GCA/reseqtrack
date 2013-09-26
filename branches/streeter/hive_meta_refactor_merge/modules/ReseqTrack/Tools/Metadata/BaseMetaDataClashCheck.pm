package ReseqTrack::Tools::Metadata::BaseMetaDataClashCheck;

use strict;
use warnings;
use base qw(ReseqTrack::Tools::Metadata::BaseMetadataAddIn);

sub report {
	my ($self) = @_;
	my $dbc = $self->reseq_db()->dbc();

	my $attr_clashes_count = $self->report_attribute_clashes();
	my $attr_col_clashes_count = $self->report_attribute_column_clashes();

}

sub report_attribute_clashes {
	my ($self) = @_;
	my $dbc    = $self->reseq_db()->dbc();
	my $sth    = $dbc->prepare( $self->_sql_attribute_2_attribute_clash );

	$sth->execute();

	my $white_list     = $self->attribute_2_attribute_whitelist;
	my $printed_header = undef;
	my $clashes_count = 0;
	while ( my $rs = $sth->fetchrow_arrayref ) {
		my ( $attribute_name, @table_names ) = @$rs;

		my $is_clash =
			( !exists $white_list->{$attribute_name}->{ $table_names[0] }
				->{ $table_names[1] }
				&& !
				exists $white_list->{$attribute_name}->{ $table_names[1] }
				->{ $table_names[0] } );

		if ( $is_clash && !$printed_header ) {
			print 'Potential attribute naming clash' . $/;
			print join( "\t", "Attribute name", "Table name 1", "Table name 2" ) . $/;
			$printed_header = 1;
		}
		if ($is_clash) {
			print join( "\t", $attribute_name, @table_names ) . $/;
			$clashes_count++;
		}
	}
	return $clashes_count;
}

sub report_attribute_column_clashes {
	my ($self) = @_;
	my $dbc = $self->reseq_db()->dbc();

	my $dbname = $dbc->dbname();
	throw("No dbname available for information schema query") unless $dbname;

	my $sth = $dbc->prepare( $self->_sql_attribute_2_column_clash );
	$sth->execute($dbname);

	my $white_list     = $self->attribute_2_column_whitelist;
	my $printed_header = undef;
	my $clashes_count;
	
	while ( my $rs = $sth->fetchrow_arrayref ) {
		my ( $column_name, $col_table_name, $attribute_name, $attr_table_name )
			= @$rs;

		my $is_clash =
			( !exists $white_list->{$col_table_name}->{$column_name}->{$attr_table_name}->{$attribute_name} );

		if ( $is_clash && !$printed_header ) {
			print 'Potential attribute/column naming clash' . $/;
			print join( "\t", "Table name (column)", "Column name", "Table name (attribute)", "Attribute name" ) . $/;
			$printed_header = 1;
		}
		if ($is_clash) {
			print join( "\t", $col_table_name,$column_name,$attr_table_name,$attribute_name ) . $/;
			$clashes_count++;
		}
	}
	return $clashes_count;
}

sub attribute_2_attribute_whitelist {
	return {
		'attribute_name' => {
			'table_name_1' => {
				'table_name_2' => 1
			}
		}
	};
}

sub attribute_2_column_whitelist {
	return {
		'col_table_name' => {
			'column_name' => {
				'attribute_table_name' => {
					'attribute_name' => 1 
				}
			}
		}
	};		
}

sub _sql_attribute_2_attribute_clash {
	return <<'END_SQL';
select a.attribute_name, least(a.table_name,b.table_name),greatest(a.table_name,b.table_name) from 
(select distinct table_name, attribute_name from attribute) a,
(select distinct table_name, attribute_name from attribute) b
where a.table_name <> b.table_name
  and a.attribute_name = b.attribute_name 
END_SQL
}

sub _sql_attribute_2_column_clash {
	return <<'END_SQL';
select c.column_name, c.table_name, a.attribute_name, a.table_name from
information_schema.columns c, 
(select distinct table_name, attribute_name from attribute) a
where c.table_schema = ?
  and c.table_name in ('run','sample','experiment','study')
  and a.table_name in ('run','sample','experiment','study')
  and a.attribute_name = c.column_name;
END_SQL
}
1;