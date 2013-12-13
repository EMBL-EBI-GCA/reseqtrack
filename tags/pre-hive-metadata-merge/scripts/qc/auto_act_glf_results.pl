#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::CollectionUtils;
use ReseqTrack::DBSQL::GenotypeResultsAdaptor;
use ReseqTrack::Tools::FileSystemUtils
  qw (create_tmp_process_dir delete_directory get_lines_from_file );
use ReseqTrack::Tools::QC::GLFTools;
use ReseqTrack::Tools::GeneralUtils qw (get_params get_open_file_handle);
use Getopt::Long;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use File::Basename;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::VerifyBamIDReadGroup;
use ReseqTrack::Tools::RunVerifyBamID;
use ReseqTrack::Tools::ArchiveUtils qw(move_archive_files_with_move_list);
my %input;
my @match_verify_bam_swap;
my @new_sample_designation;
my %processed_run_id;    # track what run_ids have been processed already

my $ftp       = '/nfs/1000g-archive/vol1/ftp/';
my $withdrawn = '/nfs/1000g-archive/vol1/withdrawn/';

$input{run} = 0;

GetOptions(
	\%input,                    'dbhost=s',
	'dbname=s',                 'dbuser=s',
	'dbpass=s',                 'dbport=s',
	'working_dir=s',            'verbose!',
	'skip_platform=s',          'debug!',
	'phase1_suppressed_file=s', 'seq_index=s',
	'force_withdraw=s',         'force_reinstate=s',
	'force_keep=s',             'run!',
);

$input{phase1_suppressed_file} = './20110304.run_level_qc.csv';

if ( defined $input{cfg_file} ) {
	get_params( $input{cfg_file}, \%input );
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
	-host   => $input{dbhost},
	-user   => $input{dbuser},
	-port   => $input{dbport},
	-dbname => $input{dbname},
	-pass   => $input{dbpass},
);
my $force_withdraw_runs = get_lines_from_file( $input{force_withdraw} )
  if ( defined $input{force_withdraw} );
my $force_reinstate_runs = get_lines_from_file( $input{force_reinstate} )
  if ( defined $input{force_reinstate} );
my $keep_in_place_runs = get_lines_from_file( $input{force_keep} )
  if ( defined $input{force_keep} );

# get read group flagged as failed VerifyBamID form DB
my $db_failed_status = get_failed_run_ids_from_VerifyBamIDReadGroup($db);

if ( defined $input{seq_index} ) {
	$db_failed_status =
	  get_failed_run_ids_from_seq_index( $input{seq_index}, $db_failed_status );
}

# read PHASE 1 correction files bases on VerifyBamID results
my $phase1_hash = inhale_phase_1_file( $input{phase1_suppressed_file} );

#Check run_meta_info & genotype_results tables to see if
# any sample swaps have been corrected by sequencing centers.
my $table_swaps = look_for_sample_swaps_in_tables($db);

#Check if sample swaps that have been corrected match VerifyBamID results.
( @match_verify_bam_swap, @new_sample_designation ) =
  check_sample_swaps_match_verifybam( $phase1_hash, $table_swaps );

#Go through genotype_results tables and get collection_ids of failed runs
my $failed_gtc_run_ids = get_genotype_check_run_ids_by_verdict( $db, "FAILED" );

#Go through genotype_results tables and get collection_ids of passed runs
my $passed_gtc_run_ids = get_genotype_check_run_ids_by_verdict( $db, "PASSED" );

#Get list of runs that are active but have "FAILED" entry in
#genotype_results table.

print "\nGet list of runs that are active but have 'FAILED' entry in\n";
print "genotype_results table.\n";

my $collections_to_withdraw =
  get_collections_to_process( $db, $failed_gtc_run_ids, "FILTERED_FASTQ" );

#Get list of runs that are withdrawn but have now have "PASSED" entry in
#genotype_results table.

print "\nGet list of runs that are withdrawn but have now have 'PASSED'\n";
print "entry in genotype_results table\n";

my $collections_to_reinstate =
  get_collections_to_process( $db, $passed_gtc_run_ids,
	"WITHDRAWN_FILTERED_FASTQ" );

#do not touch Phase 1 VerifyBamId baddies. Remove them for processing hashes

print "\nProcessing against VerifyBamID results in 20110304.run_level_qc.csv\n";

remove_phase_1_SUPPRESSED_from_coll_hashes( $collections_to_withdraw,
	$phase1_hash, " flagged to withdraw" );

remove_phase_1_SUPPRESSED_from_coll_hashes( $collections_to_reinstate,
	$phase1_hash, " flagged to reinstate" );

# $Now check what is in DB VerifyBamIDReadGroup table. If bad leave alone

print "\nProcessing against VerifyBamID results in DB\n";

remove_phase_1_SUPPRESSED_from_coll_hashes( $collections_to_withdraw,
	$db_failed_status, " flagged to withdraw" );

remove_phase_1_SUPPRESSED_from_coll_hashes( $collections_to_reinstate,
	$db_failed_status, " flagged to reinstate" );

###########################################################################################
print "\n\nSTARTING FORCE SECTION\n\n";
apply_force_actions_to_collections( $db, $collections_to_withdraw,
	$force_withdraw_runs, "FILTERED_FASTQ", "withdraw" );
apply_force_actions_to_collections( $db, $collections_to_reinstate,
	$force_reinstate_runs, "WITHDRAWN_FILTERED_FASTQ", "reinstate" );

apply_keep_actions_to_collections( $collections_to_withdraw,
	$keep_in_place_runs, "withdraw" );

apply_keep_actions_to_collections( $collections_to_reinstate,
	$keep_in_place_runs, "reinstate" );
###########################################################################################

# OK should have eveything split up. Now do something with it.

#1. Withdraw run_id's with type = "FILTERED_FASTQ" and "FAILED" glf check.
#   These should be 'public' and not SUPRRESSED' in
#   20110304.run_level_qc.csv.gz file.

move_retype_runs( $collections_to_withdraw, $db, $ftp, $withdrawn,
	"WITHDRAWN_FILTERED_FASTQ", $input{run} );

move_retype_runs( $collections_to_reinstate, $db, $withdrawn, $ftp,
	"FILTERED_FASTQ", $input{run} );

#Cleanup event_complete table for moved fastq files
#select event_id  from event where name like "genotype_%withdrawn%" and type = "WITHDRAWN_FILTERED_FASTQ";
#select event_id  from event where name like "genotype_%active%" and type = "FILTERED_FASTQ";

#select * from event_complete where table_name = "collection"
# and other_id in
#( select collection_id from collection where collection.name = "SRR037772" and
# collection.type = "WITHDRAWN_FILTERED_FASTQ") and
# event_id in ( '24','25','26','27') ;
#

sub move_retype_runs {

	my ( $colls, $db, $cur_dir, $new_dir, $new_type, $run ) = @_;
	my $ctr = 0;

	print "Moving public\/non-SUPPRESSED runs to $new_dir\n";

	if ( !keys %$colls ) {
		print "No runs to withdraw\n";
		return;
	}

	my $ca = $db->get_CollectionAdaptor;
	my $fa = $db->get_FileAdaptor;
	my %move_hash;

  COLLECTION:
	foreach my $key ( keys %$colls ) {
		$ctr++;

		my $collection = $$colls{$key};

		#print "Moving files associated with $key type = ",
		#      $collection->type, "\n";

		my $fastq = $collection->others;

		#do a check if clobbering something
		foreach my $file (@$fastq) {
			my $new_loc = $file->name;
			$new_loc =~ s/$cur_dir/$new_dir/;

			if ( -e $new_loc ) {
				print
				  "$key $new_loc exists already. Skipping this collection\n";
				next COLLECTION;

			}
		}

		print "$ctr: ", scalar(@$fastq), " files to move for $key";
		print " old type: ", $collection->type;
		print " new type: $new_type\n";

		foreach my $file (@$fastq) {
			my $new_loc = $file->name;
			$new_loc =~ s/$cur_dir/$new_dir/;
			my $new_dir = dirname($new_loc);

			my $new_file = copy_file_object($file);
			$new_file->dbID( $file->dbID );
			$new_file->type($new_type);

			my $history = create_history( $new_file, $file );
			if ($history) {
				$new_file->history($history);
				$fa->update($new_file) if $run;
			}
			else {
				print STDERR
"There appears to be no difference between the two file objects\n";
			}

			$move_hash{ $file->name } = $new_loc;

			#print $file->name, " to\n$new_loc\n\n";

		}

		print Dumper %move_hash;

		#retype collection;
		my $new_col = copy_collection_object($collection);

		$new_col->type($new_type);
		my $history = create_collection_history_object( $new_col, $collection );

		if ($history) {
			$new_col->history($history);
			$ca->update_type($new_col) if $run;
		}

	}

	move_archive_files_with_move_list( $db, \%move_hash, 0, $run);
	
	
	print '=' x 50;
	print "\n";
	return;

}

sub apply_keep_actions_to_collections {
	my ( $coll_hash, $run_ids, $action ) = @_;
	my @delete_from_coll_hash;

	if ( !defined($run_ids) ) {
		print "No force action for :$action\n";
		return;
	}

	print "Have " . scalar(@$run_ids) . " to make sure are not moved\n";

	foreach my $ri (@$run_ids) {

		if ( defined $$coll_hash{$ri} ) {
			print "Removing $ri from $action\n";
			push( @delete_from_coll_hash, $ri );
		}
	}

	my $total = scalar(@delete_from_coll_hash);

	print "Total to remove for collhash = $total\n";
	print "===============================================\n";

	foreach my $del_run_id (@delete_from_coll_hash) {
		remove( $$coll_hash{$del_run_id} );
	}

	return $coll_hash;
}

#================================================
sub apply_force_actions_to_collections {
	my ( $db, $coll_hash, $run_ids, $type, $action ) = @_;

	if ( !defined($run_ids) ) {
		print "No force action for :$action\n";
		return;
	}
	my $ca = $db->get_CollectionAdaptor;

	print "Have " . scalar(@$run_ids) . " to check\n";

	foreach my $ri (@$run_ids) {

		if ( defined $$coll_hash{$ri} ) {
			print "$ri already flagged for $action\n";
			next;
		}

		my $collection = $ca->fetch_by_name_and_type( $ri, "$type" );

		if ($collection) {
			print "Adding $ri ($type) to force $action set\n";
			$$coll_hash{$ri} = $collection;
		}
		else {
			print "$ri $type does not exist\n";
		}
	}

	my $total = scalar( keys %$coll_hash );

	print "Total to $action = $total\n";
	print "===============================================\n";

	#  print Dumper %$coll_hash;

	return $coll_hash;
}

sub get_failed_run_ids_from_seq_index {
	my ( $file, $tmp ) = @_;

	open my $IN, '<', $file || "No read on $file\n";
	print "Loading from $file\n";
	my $ctr = 0;
	while (<$IN>) {

		chomp;

		my @aa = split /\t/;

		if ( $aa[22] =~ /verify/i || $aa[22] =~ /related/i ) {

			if ( defined $$tmp{ $aa[2] } ) {

				#	print "Already got $aa[2] ($ctr)\n";
				next;
			}

			$ctr++;
			$$tmp{ $aa[2] }{run_status} = "SUPPRESS";

		}
	}
	close $IN;

	my @bad_run_ids = keys %$tmp;

	print "Now have verify/related run_id entries in hash "
	  . scalar @bad_run_ids
	  . " to update\n";
	print "Added $ctr failed VerifyBamID run_ids from index file\n\n";

	sleep 3;

	return $tmp;
}

sub get_failed_run_ids_from_VerifyBamIDReadGroup {
	my $db     = shift;
	my $failed = 1;

	my $vra = $db->get_VerifyBamIDReadGroupAdaptor();

	my ( $objs, $run_ids ) = $vra->fetch_by_status($failed);

	print "Got ", scalar( keys(%$run_ids) ), " run_ids flagged as '1'";
	print " in VerifyBamIDReadGroup\n\n";

	foreach my $key ( keys %$run_ids ) {

		#    print $key,"\n";
		$$run_ids{$key}{run_status} = "SUPPRESS";
	}

	return \%$run_ids;

}

sub check_sample_swaps_match_verifybam {
	my ( $phase1_hash, $table_swaps ) = @_;
	my $ctr = 0;

	my @matches_verifybam_swap;
	my @new_sample_designation;

	foreach my $key ( keys %$table_swaps ) {

		# print $key,"\n";
		if ( defined $$phase1_hash{$key} ) {

			my $status = $$phase1_hash{$key}{sample_status};

			#next if ( $status ne "SAMPLE_SWAP" );

			if ( defined $$table_swaps{$key}{rmi_sample_run} ) {
				if ( $$phase1_hash{$key}{reason} eq
					$$table_swaps{$key}{rmi_sample_run} )
				{
					print "Corrected sample swap matches VerBamID $key was ";
					print $$phase1_hash{$key}{sample},         " now ";
					print $$phase1_hash{$key}{reason},         "  ";
					print $$table_swaps{$key}{rmi_sample_run}, "\n";
					push( @matches_verifybam_swap, $key );
					next;
				}
			}
		}

		push( @new_sample_designation, $key );

		#print "new sample id for $key:";
		#print " gt table = ",$$table_swaps{$key}{gt_claimed}, " ";
		#print " rmi table= ",$$table_swaps{$key}{rmi_sample_run} ,"\n";

	}

	print scalar(@matches_verifybam_swap),
	  " run_ids match VerifyBamID sample swap\n";
	print scalar(@new_sample_designation),
	  " run_ids have new sample designation\n";

	return ( @matches_verifybam_swap, @new_sample_designation );
}

sub remove_phase_1_SUPPRESSED_from_coll_hashes {

	my ( $coll_hash, $p1_hash, $comment ) = @_;
	my $ctr           = 0;
	my $failed_locked = 0;
	my @to_remove_from_processing;

	my $total = scalar( keys %$coll_hash );
	print "\nStarting total $comment $total ::  \n";

	foreach my $ri ( keys %$coll_hash ) {
		$ctr++;

		print "$ctr $ri: ";
		if ( defined $$p1_hash{$ri}{run_status} ) {

			#print Dumper $$p1_hash{$ri};
			if (   $$p1_hash{$ri}{run_status} eq "SUPPRESS"
				|| $$p1_hash{$ri}{run_status} eq "EXCLUDE" )
			{
				print " failed QC process";

				#	print "P1 $ri :", $$p1_hash{$ri}{run_status},"\n";
				$failed_locked++;
				push( @to_remove_from_processing, $ri );
			}
		}
		print "\n";
	}
	print "Failed VerBamID = $failed_locked  :: removing ",
	  scalar @to_remove_from_processing, "  ";

	foreach my $drop (@to_remove_from_processing) {
		delete( $$coll_hash{$drop} );
	}

	$total = scalar( keys %$coll_hash );
	print "run_ids to process = $total\n";

}

sub get_collections_to_process {
	my ( $db, $failed_run_ids, $type ) = @_;
	my %coll_hash;

	my $ca = $db->get_CollectionAdaptor;

	foreach my $ri (@$failed_run_ids) {

		my $collection = $ca->fetch_by_name_and_type( $ri, "$type" );
		if ($collection) {
			$coll_hash{$ri} = $collection;
		}
	}

	my $total = scalar( keys %coll_hash );

	print "Found $total collections of type $type to process\n\n";

	return \%coll_hash;
}

sub get_genotype_check_run_ids_by_verdict {
	my ( $db, $verdict ) = @_;

	my $sql =
	    "select gt.name, gt.verdict, rmi.run_id,rmi.status "
	  . " from  genotype_results gt , run_meta_info rmi  where "
	  . "gt.verdict = \"$verdict\" and gt.name = rmi.run_id and "
	  . " rmi.status = \"public\" ";

	my @results;

	my $sth = $db->dbc->prepare($sql);
	$sth->execute();
	my $ctr = 0;
	while ( my $rowArrayref = $sth->fetchrow_arrayref ) {
		$ctr++;

		#print "$ctr @$rowArrayref[0] @$rowArrayref[1]\n";
		push( @results, @$rowArrayref[0] );
	}
	print "Got ", scalar(@results), " runs that $verdict genotype check\n";

	return \@results;
}

sub inhale_phase_1_file {

	#Do not reinstate anything from this file, unless sample swap has been
	#corrected.

	my $p1_file = shift;

	my %p1_hash;
	die "++No file found +$p1_file+" if ( !-e $p1_file );

	my $FH = get_open_file_handle($p1_file);

	print "Reading $p1_file\n";
	while (<$FH>) {
		next if (/^run_id/);

		chomp;
		my (
			$run_id,     $center,   $sample,
			$population, $platform, $has_phase1_bam,
			$run_status, $reason,   $sample_status
		) = split /\,/;

		$p1_hash{$run_id}{run_id}         = $run_id;
		$p1_hash{$run_id}{center}         = $center;
		$p1_hash{$run_id}{sample}         = $sample;
		$p1_hash{$run_id}{population}     = $population;
		$p1_hash{$run_id}{platform}       = $platform;
		$p1_hash{$run_id}{has_phase1_bam} = $has_phase1_bam;
		$p1_hash{$run_id}{run_status}     = $run_status;
		$p1_hash{$run_id}{reason}         = $reason;
		$p1_hash{$run_id}{sample_status}  = $sample_status;
		if ( $run_status ne "OK" ) {

# print "$run_id ",$p1_hash{$run_id}{run_status}, " \t", $p1_hash{$run_id}{reason}, "\n";
		}
	}
	my $swap_total = 0;
	foreach my $i ( keys %p1_hash ) {
		if ( $p1_hash{$i}{reason} =~ /TO\:/ ) {
			$swap_total++;
			$p1_hash{$i}{reason} =~ s/^TO\://;

	   #      print $p1_hash{$i}{run_status}, " \t", $p1_hash{$i}{reason}, "\n";
		}

	}
	print "Found $swap_total run_ids with sample swaps in $p1_file\n";

	return \%p1_hash;
}

sub look_for_sample_swaps_in_tables {
	my $db = shift;
	my @results;
	my %table_swaps;

	my $sql = "select gt.name, gt.claimed, rmi.run_id, rmi.sample_name from
  genotype_results gt , run_meta_info rmi 
  where gt.name = rmi.run_id and gt.claimed != rmi.sample_name
  and  rmi.run_id in (select run_id from run_meta_info where status = \"public\")";

	#print $sql,"\n";

	my $sth = $db->dbc->prepare($sql);
	$sth->execute();

	while ( my $rowArrayref = $sth->fetchrow_arrayref ) {

		#print "@$rowArrayref\n";
		push( @results, $rowArrayref );
		$table_swaps{ @$rowArrayref[0] }{gt_claimed}     = @$rowArrayref[1];
		$table_swaps{ @$rowArrayref[0] }{rmi_sample_run} = @$rowArrayref[3];

		# print " gt @$rowArrayref[1] rmi @$rowArrayref[3]\n";
	}

	print "Got ", scalar @results;
	print
	  " cases of sample swap(s) in run_meta_info\/genotype_results tables\n";

	#print Dumper %table_swaps;
	#exit;
	return \%table_swaps;
}

