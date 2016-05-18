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
use ReseqTrack::Tools::AttributeUtils
  qw(create_attribute_for_object remove_outdated_attributes);

use feature 'switch';

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

throw(
"No work to do: specify -load_new and/or -update_existing. Use -help to see options"
) unless ( $load_new || $update_existing );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(%db_params);
$db->dbc->db_handle->{'AutoCommit'} = 0;

# the order of types matters - e.g. experiments reference studies, runs reference samples and experiments
my @types = qw(studies samples experiments runs);

TYPE: for my $type (@types) {
  my $file_path = $files{$type};

  next TYPE unless $file_path;

  print "Reading $type file: $file_path $/" if ($verbose);
  my $objects = load_from_file( $file_path, $type, $separator, $verbose );

  next TYPE unless $objects;
  print "Loading/updating $type$/" if ($verbose);

  my $adaptor          = adaptors( $db, $type );
  my $fk_handler_sub   = fk_handlers($type);
  my $allow_attributes = have_attributes($type);

  for my $object (@$objects) {
    my $current_record = $adaptor->fetch_by_source_id( $object->source_id );

    $fk_handler_sub->( $object, $db );

    if ( $update_existing && defined $current_record ) {
      print "Updating $type " . $object->source_id() . $/ if ($verbose);
      $object->dbID( $current_record->dbID() );
      $adaptor->update($object);

      if ($allow_attributes) {
        $adaptor->store_attributes( $object, 1 );
        remove_outdated_attributes( $object, $current_record,
          $db->get_AttributeAdaptor );
      }

    }
    if ( $load_new && !defined $current_record ) {
      print "Storing $type " . $object->source_id() . $/ if ($verbose);
      $adaptor->store($object);
      $adaptor->store_attributes($object) if ($allow_attributes);
    }
  }

}

$db->dbc->db_handle->commit();

sub load_from_file {
  my ( $file, $type, $separator, $verbose, ) = @_;

  my $columns          = columns($type);
  my $allow_attributes = have_attributes($type);
  my $class            = classes($type);

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
            create_attribute_for_object( $object, $attribute_name,
            $attribute_value );
        }
        $object->attributes( \@attr );
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

  my $count = scalar @objects;
  print "Found $count entries $/" if ($verbose);

  return \@objects;
}

sub print_template {
  my ( $separator, $files ) = @_;

  while ( my ( $file_type, $file_path ) = each %$files ) {
    next unless $file_path;

    throw("Template path for $file_type already exists: $file_path ")
      if ( -e $file_path );

    my @columns = @{ columns($file_type) };
    throw("Could not find columns for $file_type template") unless @columns;

    push @columns, "#add attribute columns here"
      if ( have_attributes($file_type) );

    open( my $fh, '>', $file_path )
      or throw("Could not write to $file_path: $!");

    print $fh join( $separator, @columns ) . $/;

    close $fh;
  }
}

# columns to be used in the input files. used both for parsing and for generating templates
#first column in these lists must always map to the source_id column
sub columns {
  my ($t) = @_;
  given ($t) {
    when ('studies') { return [qw(STUDY_ID	TITLE	TYPE	STATUS)] }
    when ('samples') {
      return [
        qw(SAMPLE_ID	SAMPLE_TITLE	STATUS	CENTER_NAME	SAMPLE_ALIAS	TAX_ID	SCIENTIFIC_NAME	COMMON_NAME	ANONYMIZED_NAME	INDIVIDUAL_NAME	)
        ]
    }
    when ('experiments') {
      return [
        qw(EXPERIMENT_ID	STUDY_ID	SAMPLE_ID	STATUS	CENTER_NAME	EXPERIMENT_ALIAS	INSTRUMENT_PLATFORM	INSTRUMENT_MODEL	LIBRARY_LAYOUT	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	PAIRED_NOMINAL_LENGTH	PAIRED_NOMINAL_SDEV)
        ]
    }
    when ('runs') {
      return [
        qw(RUN_ID	EXPERIMENT_ID	RUN_ALIAS	STATUS	CENTER_NAME	RUN_CENTER_NAME	INSTRUMENT_PLATFORM	INSTRUMENT_MODEL)
        ]
    }
  }
}

# sub-routines to turn the supplied ids for related entries in to db IDs. Also responsible for throwing errors if the values aren't set or don't relate to anything in the db
sub fk_handlers {
  my ($t) = @_;
  given ($t) {
    when ('studies') {
      return sub { }
    }
    when ('samples') {
      return sub { }
    }
    when ('experiments') {
      return sub {
        my ( $e, $db ) = @_;

        throw( "No study ID set for experiment " . $e->source_id() )
          unless ( $e->study_id );
        throw( "No sample ID set for run " . $r->source_id() )
          unless ( $r->sample_id );

        my $st_a = $db->get_StudyAdaptor();
        my $stid = $e->study_id();
        my $st   = $st_a->fetch_by_source_id($stid);
        throw( "Could not find study $stid for experiment " . $e->source_id() )
          unless $st;
        $e->study_id( $st->dbID() );

        my $sa_a = $db->get_SampleAdaptor();
        my $said = $e->sample_id();
        my $sa   = $sa_a->fetch_by_source_id($said);
        throw( "Could not find sample $said for experiment " . $e->source_id() )
          unless $sa;
        $e->sample_id( $sa->dbID() );

        }
    }
    when ('runs') {
      return sub {
        my ( $r, $db ) = @_;

        throw( "No experiment ID set for run " . $r->source_id() )
          unless ( $r->experiment_id );

        my $ea  = $db->get_ExperimentAdaptor();
        my $eid = $r->experiment_id();
        my $e   = $ea->fetch_by_source_id($eid);
        throw( "Could not find experiment $eid for run " . $r->source_id() )
          unless $e;
        $r->experiment_id( $e->dbID() );
        }
    }
  }
}

# which data types can use the attribute table each data type
sub have_attributes {
  my ($t) = @_;
  given ($t) {
    when ('studies')     { return 0 }
    when ('samples')     { return 1 }
    when ('experiments') { return 1 }
    when ('runs')        { return 0 }
  }
}

# the class to use when creating each data type
sub classes {
  my ($t) = @_;
  given ($t) {
    when ('studies')     { return 'ReseqTrack::Study' }
    when ('samples')     { return 'ReseqTrack::Sample' }
    when ('experiments') { return 'ReseqTrack::Experiment' }
    when ('runs')        { return 'ReseqTrack::Run' }
  }
}

# the db adaptors required for each data type
sub adaptors {
  my ( $db, $t ) = @_;
  given ($t) {
    when ('studies')     { return $db->get_StudyAdaptor() }
    when ('samples')     { return $db->get_SampleAdaptor() }
    when ('experiments') { return $db->get_ExperimentAdaptor() }
    when ('runs')        { return $db->get_RunAdaptor() }
  }
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
        -verbose, print excessive but not quite useful information to STDOUT while running

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
