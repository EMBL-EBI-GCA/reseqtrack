use strict;
use warnings;
use Data::Dumper;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Study;
use ReseqTrack::Sample;
use ReseqTrack::Experiment;
use ReseqTrack::Run;

use Getopt::Long;

my %db_params;
my %files;
my $write_templates;
my $load_new;
my $update_existing;
my $help;
my $separator = "\t";

&GetOptions(
	'dbhost=s'         => \$db_params{-host},
	'dbname=s'         => \$db_params{-dbname},
	'dbuser=s'         => \$db_params{-user},
	'dbpass=s'         => \$db_params{-pass},
	'dbport=s'         => \$db_params{-port},
	'studies=s'        => \$files{studies},
	'samples=s'        => \$files{samples},
	'experiments=s'    => \$files{experiments},
	'runs=s'           => \$files{runs},
	'write_templates!' => \$write_templates,
	'load_new!'        => \$load_new,
	'update_existing!' => \$update_existing,
	'help!'            => \$help,
	'separator=s'      => \$separator,
);

throw("Cannot write templates and read files")
	if ( $write_templates && ( $load_new || $update_existing ) );

if ($help) {
	exec( 'perldoc', $0 );
	exit(0);
}
if ($write_templates) {
	print_template( $separator, \%files );
	exit;
}

throw("No work to do: specify -load_new and/or -update_existing")
	unless ( $load_new || $update_existing );

my $meta_data_in_file = load_meta_data_from_files( \%files, $separator );

print Dumper($meta_data_in_file);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(%db_params);
$db->dbc->db_handle->{'AutoCommit'} = 0;

my $adaptors = adaptors($db);

for my $type ( keys %$meta_data_in_file ) {
	my $adaptor = $adaptors->{$type};
	my $objects = $meta_data_in_file->{$type};
	
	for my $object (@$objects) {
		my $current_record = $adaptor->fetch_by_dbID( $object->dbID );
	
		if ($update_existing && defined $current_record  ){
			$adaptor->update($object);
		}	
		if ($load_new && !defined $current_record ){
			$adaptor->store($object);
		}
	}
	
}

$db->dbc->db_handle->commit();

sub load_meta_data_from_files {
	my ( $files, $separator ) = @_;

	my $columns        = columns();
	my $has_attributes = have_attributes();
	my $classes        = classes();

	my %meta_data_in_file;
	while ( my ( $file_type, $file_path ) = each %$files ) {
		next unless $file_path;

		my $cols             = $columns->{$file_type};
		my $allow_attributes = $has_attributes->{$file_type};
		my $class            = $classes->{$file_type};

		$meta_data_in_file{$file_type} =
			load_from_file( $file_path, $cols, $allow_attributes, $class,
			$separator );
	}

	return \%meta_data_in_file;
}

sub load_from_file {
	my ( $file, $columns, $allow_attributes, $class, $separator ) = @_;

	my @objects;
	my @header;
	my %required_cols;
	map { $required_cols{$_} = 1 } @$columns;

	open( my $fh, '<', $file )
		or throw("Could not read from $file: $!");

	while (<$fh>) {
		chomp;
		my @vals = split $separator;

		if (@header) {
			my %item;
			my %attributes;

			for ( my $i = 0 ; $i < scalar(@header) ; $i++ ) {
				my $h = $header[$i];
				my $v = $vals[$i];

				if ($h eq $columns->[0] ){
					$item{-dbID} = $v;
				}
				elsif ( $required_cols{$h} ) {
					$item{ '-' . $h } = $v;
				}
				elsif ($allow_attributes) {
					$attributes{$h} = $v;
				}
				else {
					throw(
"Unexpected outcome: $h is not a standard columnn in $file and this object does not allow arbitrary objects"
					);
				}
			}

			my $object = $class->new(%item);
			if ($allow_attributes) {
				$object->attributes( \%attributes );
			}
			push @objects, $object;
		}
		else {
			@header = @vals;

			# check that the headers are valid for this file
			unless ($allow_attributes) {
				for my $h (@header) {
					throw( "Column $h is not an expected column in $file. We expect "
							. join( "\t", @$columns ) )
						unless $required_cols{$h};
				}
			}
		}

	}

	close $fh;

	return \@objects;
}

sub print_template {
	my ( $separator, $files ) = @_;

	my $file_columns    = columns();
	my $have_attributes = have_attributes();

	while ( my ( $file_type, $file_path ) = each %$files ) {
		next unless $file_path;

		throw("Template path for $file_type already exists: $file_path ")
			if ( -e $file_path );

		my @columns = @{ $file_columns->{$file_type} };
		push @columns, "#add attribute columns here"
			if $have_attributes->{$file_type};

		throw( "Could not find columns for $file_type template. Lookup table: "
				. Dumper($file_columns) )
			unless @columns;

		open( my $fh, '>', $file_path )
			or throw("Could not write to $file_path: $!");

		print $fh join( $separator, @columns ) . $/;

		close $fh;
	}
}

sub columns {
	#first column in these lists must always be the dbID
	return {
		studies => [qw(STUDY_ID	TITLE	TYPE	STATUS)],
		samples => [
			qw(SAMPLE_ID	SAMPLE_TITLE	STATUS	CENTER_NAME	SAMPLE_ALIAS	TAX_ID	SCIENTIFIC_NAME	COMMON_NAME	ANONYMIZED_NAME	INDIVIDUAL_NAME	)
		],
		experiments => [
			qw(EXPERIMENT_ID	STUDY_ID	SAMPLE_ID	STATUS	CENTER_NAME	EXPERIMENT_ALIAS	INSTRUMENT_PLATFORM	INSTRUMENT_MODEL	LIBRARY_LAYOUT	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	PAIRED_NOMINAL_LENGTH	PAIRED_NOMINAL_SDEV)
		],
		runs => [
			qw(RUN_ID	EXPERIMENT_ID	RUN_ALIAS	STATUS	CENTER_NAME	RUN_CENTER_NAME	INSTRUMENT_PLATFORM	INSTRUMENT_MODEL)
		],
	};
}

sub have_attributes {
	return {
		studies     => 0,
		samples     => 1,
		experiments => 1,
		runs        => 0
	};
}

sub classes {
	return {
		studies     => 'ReseqTrack::Study',
		samples     => 'ReseqTrack::Sample',
		experiments => 'ReseqTrack::Experiment',
		runs        => 'ReseqTrack::Run',
	};
}

sub adaptors {
	my ($db) = @_;
	return {
	studies     => $db->get_StudyAdaptor(),
	samples     => $db->get_SampleAdaptor(),
	experiments => $db->get_ExperimentAdaptor(),
	runs        => $db->get_RunAdaptor(),
	};
}

=pod

=head1 NAME

ReseqTrack/scripts/metadata/load_from_file.pl

=head1 SYNOPSIS



=head2 OPTIONS

      database options:

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on
       

      other options:

    

        -help, flag to print this help and exit


=head1 Example:


   

=cut
