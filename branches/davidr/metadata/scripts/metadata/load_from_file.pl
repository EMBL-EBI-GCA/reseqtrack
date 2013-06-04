#! /usr/bin/env perl
use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;

use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Study;
use ReseqTrack::Sample;
use ReseqTrack::Experiment;
use ReseqTrack::Run;
use ReseqTrack::Tools::StatisticsUtils qw(create_statistic_for_object);

my %db_params;
my %files;
my $write_templates;
my $load_new;
my $update_existing;
my $help;
my $separator = "\t";
my $verbose;

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
	'verbose!'         => \$verbose,
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

throw("No work to do: specify -load_new and/or -update_existing. Use -help to see options")
	unless ( $load_new || $update_existing );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(%db_params);
$db->dbc->db_handle->{'AutoCommit'} = 0;

my $meta_data_in_file =
	load_meta_data_from_files( \%files, $separator, $verbose );

my $adaptors          = adaptors($db);
my $attribute_adaptor = $db->get_AttributeAdaptor;
my $fk_handlers       = fk_handlers();
my $have_attributes   = have_attributes();

my @types = qw(studies samples experiments runs)
	; # the order of types matters - e.g. experiments reference studies, runs reference samples and experiments

TYPE: for my $type (@types) {
	my $objects = $meta_data_in_file->{$type};
	next TYPE unless $objects;

	my $adaptor          = $adaptors->{$type};
	my $fk_handler       = $fk_handlers->{$type};
	my $allow_attributes = $have_attributes->{$type};

	print STDERR "Loading/updating $type$/" if ($verbose);

	for my $object (@$objects) {
		my $current_record = $adaptor->fetch_by_source_id( $object->source_id );

		$fk_handler->( $object, $db );

		if ( $update_existing && defined $current_record ) {
			print STDERR "Updating $type " . $object->source_id() . $/ if ($verbose);
			$object->dbID( $current_record->dbID() );
			$adaptor->update($object);

			if ($allow_attributes) {
				$adaptor->store_attributes( $object, 1 );
				remove_outdated_attributes( $object, $current_record,
					$attribute_adaptor, $verbose );
			}

		}
		if ( $load_new && !defined $current_record ) {
			print STDERR "Storing $type " . $object->source_id() . $/ if ($verbose);
			$adaptor->store($object);
			$adaptor->store_attributes($object) if ($allow_attributes);
		}
	}

}

$db->dbc->db_handle->commit();

sub remove_outdated_attributes {
	my ( $new_version, $old_version, $attribute_adaptor, $verbose ) =
		@_;

	my %new_attr_names;
	map { $new_attr_names{ $_->attribute_name() } = 1 }
		@{ $new_version->statistics() };

	for my $attr ( @{ $old_version->statistics() } ) {
		my $key = $attr->attribute_name();
		if ( !exists $new_attr_names{$key} ) {
			print STDERR "Attribute $key to be deleted$/" if ($verbose);
			$attribute_adaptor->delete($attr);
		}
	}

}

sub load_meta_data_from_files {
	my ( $files, $separator, $verbose ) = @_;

	my $columns        = columns();
	my $has_attributes = have_attributes();
	my $classes        = classes();

	my %meta_data_in_file;
	while ( my ( $file_type, $file_path ) = each %$files ) {
		next unless $file_path;

		my $cols             = $columns->{$file_type};
		my $allow_attributes = $has_attributes->{$file_type};
		my $class            = $classes->{$file_type};

		print STDERR "Reading $file_type file: $file_path $/" if ($verbose);

		$meta_data_in_file{$file_type} =
			load_from_file( $file_path, $cols, $allow_attributes, $class,
			$separator );

		my $count = scalar @{ $meta_data_in_file{$file_type} };
		print STDERR "Found $count entries $/" if ($verbose);
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

				if ( $h eq $columns->[0] ) {
					$item{-SOURCE_ID} = $v;
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
			if (%attributes) {
				my @attr;
				while ( my ( $attribute_name, $attribute_value ) = each %attributes ) {
					push @attr,
						create_statistic_for_object( $object, $attribute_name,
						$attribute_value );
				}
				$object->statistics( \@attr );
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

# columns to be used in the input files. used both for parsing and for generating templates
#first column in these lists must always map to the source_id column
sub columns {
	return {
		studies => [qw(STUDY_ID	TITLE	TYPE	STATUS)],
		samples => [
			qw(SAMPLE_ID	SAMPLE_TITLE	STATUS	CENTER_NAME	SAMPLE_ALIAS	TAX_ID	SCIENTIFIC_NAME	COMMON_NAME	ANONYMIZED_NAME	INDIVIDUAL_NAME	)
		],
		experiments => [
			qw(EXPERIMENT_ID	STUDY_ID	STATUS	CENTER_NAME	EXPERIMENT_ALIAS	INSTRUMENT_PLATFORM	INSTRUMENT_MODEL	LIBRARY_LAYOUT	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	PAIRED_NOMINAL_LENGTH	PAIRED_NOMINAL_SDEV)
		],
		runs => [
			qw(RUN_ID	EXPERIMENT_ID	SAMPLE_ID	RUN_ALIAS	STATUS	CENTER_NAME	RUN_CENTER_NAME	INSTRUMENT_PLATFORM	INSTRUMENT_MODEL)
		],
	};
}

# sub-routines to turn the supplied ids for related entries in to db IDs. Also responsible for throwing errors if the values aren't set or don't relate to anything in the db
sub fk_handlers {
	return {
		studies     => sub { },
		samples     => sub { },
		experiments => sub {
			my ( $e, $db ) = @_;

			throw( "No study ID set for experiment " . $e->source_id() )
				unless ( $e->study_id );

			my $st_a = $db->get_StudyAdaptor();
			my $stid = $e->study_id();
			my $st   = $st_a->fetch_by_source_id($stid);
			throw( "Could not find study $stid for experiment " . $e->source_id() )
				unless $st;
			$e->study_id( $st->dbID() );

		},
		runs => sub {
			my ( $r, $db ) = @_;

			throw( "No experiment ID set for run " . $r->source_id() )
				unless ( $r->experiment_id );
			throw( "No sample ID set for run " . $r->source_id() )
				unless ( $r->sample_id );

			my $ea  = $db->get_ExperimentAdaptor();
			my $eid = $r->experiment_id();
			my $e   = $ea->fetch_by_source_id($eid);
			throw( "Could not find experiment $eid for run " . $r->source_id() )
				unless $e;
			$r->experiment_id( $e->dbID() );

			my $sa_a = $db->get_SampleAdaptor();
			my $said = $r->sample_id();
			my $sa   = $sa_a->fetch_by_source_id($said);
			throw( "Could not find sample $said for run " . $r->source_id() )
				unless $sa;
			$r->sample_id( $sa->dbID() );
		},
		}

}

# which data types can use the attribute table each data type
sub have_attributes {
	return {
		studies     => 0,
		samples     => 1,
		experiments => 1,
		runs        => 0
	};
}

# the class to use when creating each data type
sub classes {
	return {
		studies     => 'ReseqTrack::Study',
		samples     => 'ReseqTrack::Sample',
		experiments => 'ReseqTrack::Experiment',
		runs        => 'ReseqTrack::Run',
	};
}

# the db adaptors required for each data type
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

Load meta data into a reseqtrack database from files.

=head2 OPTIONS

      database options:

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on

      file options:
			
        Each of these will be used as a read or write target, depending on the mode, for that data type
				
        -studies
        -samples
        -experiments
        -runs

      mode options:
		
        -write_templates, will write templates for each data type to the files specified. Will not overwrite existing files.
        -load_new, read the files specified and load new entries
        -update_existing, read the files specified and update existing entries. The files will be treated as authoritative - any attributes hat are in the db, but not in the file will be removed.
				    
      misc:

        -help, flag to print this information and exit
        -separator, the column separator to use when reading or writing data files. Defaults to tab.
        -verbose, print excessive but not quite useful information to STDERR while running

=head1 Notes on usage:

All data types have standard columns. Samples and experiments can also have attributes. Any column in the meta data file that isn't in the standard list will be interpreted as an attribute.

Experiments reference studies. Runs reference samples and experiments. If the referenced entities can't be found the script will fail.

All entries are matched up comparing the source_id column in the db to the <entity>_id in the spreadsheets. 

Auto commit is disabled for this script - if an error is thrown no changes will be made.

=head1 Examples:

	Write template files for all data types:

	load_from_file.pl -write_templates -studies study_template -samples sample_template -experiments experiment_template -runs run_template
	
	Load new studies:

	load_from_file.pl -dbhost mysql-blueprint -dbport 4360 -dbuser g1krw -dbname test_metainf -dbpass mypassword  -load_new -studies study_files.tsv

	Do everything:
	
	load_from_file.pl -dbhost mysql-blueprint -dbport 4360 -dbuser g1krw -dbname test_metainf -dbpass mypassword  -load_new -update_existing -studies study_files.tsv -samples samples.tsv -experiments experiments.tsv -runs runs.tsv 
   

=cut
