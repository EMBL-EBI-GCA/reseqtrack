#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::DBSQL::DBAdaptor; 
use ReseqTrack::DBSQL::HostAdaptor;  
use ReseqTrack::Tools::BamUtils;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use File::Basename;
use File::Path;
use Getopt::Long;
use Time::Local;

$| = 1; 

my $start_run = time();

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
 	$file_type,
    $help,
    $run
);    

&GetOptions(
  'dbhost=s'     		=> \$dbhost,
  'dbname=s'     		=> \$dbname,
  'dbuser=s'     		=> \$dbuser,
  'dbpass=s'     		=> \$dbpass,
  'dbport=s'     		=> \$dbport,
  'file_type=s'			=> \$file_type,
  'help!'		 		=> \$help,
  'run!'				=> \$run,
);


my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host => $dbhost,
  -user => $dbuser,
  -port => $dbport,
  -dbname => $dbname,
  -pass => $dbpass,
    );
    
my $fa = $db->get_FileAdaptor;
my $fos = $fa->fetch_by_type($file_type);

throw("No file object is found; please check the bam file type\n") if (!$fos  || @$fos == 0);

my $ca = $db->get_CollectionAdaptor;
		
$db->dbc->disconnect_when_inactive(1);

my %collection_to_transpose;
foreach my $fo (@{$fos}) {
	my $file_path = $fo->name;			
	my ($sample, $platform, $algorithm, $project, $analysis, $chrom, $date, $pop) = CHECK_AND_PARSE_FILE_NAME($fo->name);
	next unless ($chrom =~ /mapped/i && $chrom !~ /unmapped/i);
	#next unless ($pop =~ /PUR/i ); #### FIXME, remove after testing
	my $col_name = $analysis . "_" . $pop . "_" . $platform;
	push @{$collection_to_transpose{$col_name}}, $fo;	
}

foreach my $col (keys %collection_to_transpose ) {
		
		my $type = $file_type . "_TO_TRANSPOSE";
		my $collection =  ReseqTrack::Collection->new(
		  -name => $col,
		  -others => \@{$collection_to_transpose{$col}},
		  -type => $type,  
		  -table_name => 'file',
		);
	
		if($run ) {
			$ca->store($collection); ## the store function checks for exisiting collection and update them
			print "Collection $col is loaded\n" ;
		}
		elsif (!$run) {
			print  "Collection $col is not loaded. Please set -run if you like to load them\n" ;
		}
}

#############################################################################################

=pod

=head1 NAME

create_bam_collection_to_transpose.pl

=head1 SYNOPSIS

This script organize BAM files of a given type (PHASE1_BAM, EXOME_BAM, or BAM) into collections, each collection contains BAMs from the same 
population/platform/analysis_grp. One chromosome from all BAMs in one collection would be merged into a chrN BAM for variant calling. 
	
=head1 EXAMPLE
	
1000genomes.ebi.ac.uk> perl $ZHENG_RP/bin/create_bam_collection_to_transpose.pl $WRITE_DB_ARGS -dbname zheng_var_call -file_type PHASE1_BAM -run		
