#!/usr/bin/perl -w
#script should do following in order..
# 1. Retrieve collection objects for all input IDs
# 2. Spin through all collection objects, just to check
#    all the file types are the same. Just output warning
#    if more the 1 file type associated with collection.
# 3  For each collection, change type of all files assocaited with that
#     collection ( do history objects etc). The change type of collection object.


use strict;
use warnings;

use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::CollectionUtils;
use File::Basename;
use Getopt::Long;
use Data::Dumper;

local $| = 1;

my %input;
my $all_collections;
my %run_hash;
my @name;

$input{verbose} = 0;



GetOptions(
	\%input,      'dbhost=s', 'dbname=s',   'dbuser=s',
	'dbpass=s',   'dbport=s', 'verbose!',   'list=s',
	'old_type=s', 'new_type=s', 'run','help!',
	'name=s@',
);

if($input{help}){
  useage();
}



if ( !defined $input{run} ) {
	print "\nNotice: run option off. Test mode\n\n";
	sleep(3);
}

throw "\nNo new type specified\n" if (!defined $input{new_type});
throw "\nNo old type specified\n" if (!defined $input{old_type});

my ($db,$fa,$ca) = get_db_connections(%input);

 if (defined $input{name}){
 	foreach (@{$input{name}}){
 		print "Adding $_ to process list\n";
	   $run_hash{$_} = 1;
 	}	
};

%run_hash = &get_list_of_colls($input{list} ) if ($input{list});

$all_collections = &get_collections( \%run_hash, $ca, $input{old_type} );

&coll_file_types_consistent($all_collections);


foreach my $run_id ( keys(%run_hash) ) {

	my $collection = $$all_collections{$run_id};

	print "Processing: " . $collection->name . " " . $collection->type . "\n";
	my $total = @{ $collection->others };
	my $others = $collection->others;

	print "Retyping others to $input{new_type}\n";
	foreach my $other ( @{$others} ) {
		&retype_file( $fa, $other->name, $input{new_type}, $input{run} );
	}

	my $new_col  = copy_collection_object($collection);
	my $new_type = $input{new_type};

	$new_col->type($new_type);
	my $history = create_collection_history_object( $new_col, $collection );

	if ($history) {
		$new_col->history($history);
		$ca->update_type($new_col) if $input{run};
	}
}





##### Done #####

#############################
sub get_collections {
	my ( $colls, $ca, $old_type ) = @_;

	my %collections;

	foreach my $run_id ( keys(%$colls) ) {

		my $collection = $ca->fetch_by_name_and_type( $run_id, $old_type );
		if ( !$collection ) {
			print STDERR "Failed to fetch a collection for " . $run_id
			  . " using $old_type\n";
			next;
		}

		#print "Got $run_id collection\n";
		$collections{$run_id} = $collection;
	}
	
    throw("Nothing to process") if ( !scalar( keys %collections ) );
	return \%collections;
}
#############################
sub coll_file_types_consistent {

	# check all files types are consistent, within collection.
	my ($collections) = @_;
	my $bad = 0;

	foreach my $run_id ( keys %$collections ) {

		my $this_coll = $$collections{$run_id};

		#print $run_id, "\t", $this_coll->type,"\n";

		my %type_check;
		my $others = $this_coll->others;
		foreach my $other ( @{$others} ) {

			# 	print $other->name,"\t",$other->type,  "\n";
			$type_check{ $other->type } += 1;
		}
		my $number_of_types = scalar( keys %type_check );

		#print "number of types in this collection = $number_of_types\n";

		if ( $number_of_types > 1 ) {
			warning "More than 1 type in collection = $run_id\n",;
		}
	}

	return;
}
#############################
sub get_list_of_colls {

	my $list = shift;

	if ( !-e $list ) {
		print "No file list\n";
		exit;
	}
	my $lines = get_lines_from_file($list);
	my %run_hash;

	foreach my $line (@$lines) {
		my $name = basename($line);
		$name =~ /^([E|S]RR\d+)/;
		my $run_id = $1;
		$run_hash{$run_id}  = 1 ;
	}

	print "Have ". scalar( keys %run_hash ). " collections to retype\n";
	throw("Nothing to process") if ( !scalar( keys %run_hash ) );

	return %run_hash;
}
#############################
sub retype_file {

	my ( $fa, $line, $new_type, $run ) = @_;

	my $name = basename($line);

	my $files = $fa->fetch_all_like_name($name);
	if ( !$files || @$files == 0 ) {
		print STDERR $line . " doesn't exist in the database\n";
		next;
	}

	if ( @$files >= 2 ) {
		print STDERR "Have " . @$files . " which match " . $name . "\n";
		foreach my $file (@$files) {
			print STDERR $file->name . "\t" . $file->type . ".Skipping\n";
		}
		return;
	}

	print $files->[0]->name . " " . $files->[0]->type . "\n";
	my $old_file = $files->[0];
	my $new_file = copy_file_object($old_file);
	$new_file->dbID( $old_file->dbID );
	$new_file->type($new_type);

	my $history = create_history( $new_file, $old_file );
	if ($history) {
		$new_file->history($history);
		$fa->update($new_file) if $run;
	}
	else {
		print STDERR
		  "There appears to be no difference between the two file objects\n";
	}
	return;
}



sub get_db_connections{
	my $input = shift;
	
	my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $input{dbhost},
    -user   => $input{dbuser},
    -port   => $input{dbport},
    -dbname => $input{dbname},
    -pass   => $input{dbpass},
    );

    my $fa = $db->get_FileAdaptor;
    my $ca = $db->get_CollectionAdaptor;
	
	
	
	return ($db,$fa,$ca);
	
}

sub useage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/files/retype_collection.pl

=head1 SYNOPSIS

  This script will take a list of collection ids and a specified 
  new type and alter the corresponding 'type' of all the files 
  (in 'file' table) associated with that collection object. Then
   it should change the type of the collection (collection table) 
   to the new type.

=head1 OPTIONS

  Database options

  -dbhost, the name of the mysql-host
  -dbname, the name of the mysql database
  -dbuser, the name of the mysql user
  -dbpass, the database password if appropriate
  -dbport, the port the mysql instance is running on, this defaults to 4197
     the standard port for mysql-g1kdcc.ebi.ac.uk

  Input options
  -list,     text file contain full pathnames of files that need type changing.
  -new_type, new collection type
  -old_type, current type of collection to be changed
  -name,     name of collection to be changed
  -run,      include to execute file type changes

 Example:
    perl retype_collection.pl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass XXXX
    -dbport 4197 -dbname lec_test_track -list file.list old_type BAM -new_type
    WITHDRAWN_BAM  -run
     
    perl retype_collection.pl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass XXXX
    -dbport 4197 -dbname lec_test_track -name  ERR00001 old_type FILTERED_FASTQ
    -new_type WITHDRAWN_FILTERED_FASTQ  -run
     
    

=cut

