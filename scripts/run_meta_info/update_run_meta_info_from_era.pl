#!/sw/arch/bin/perl 

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::UpdateRunMetaInfo;
use ReseqTrack::Tools::ERAUtils;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $era_dbuser;
my $era_dbpass;
my $run;
my $help;
my $store_new;
my $update_existing;
my $fix_sample;
my $withdraw_suppressed;
my $check_individual;
my $all_checks;
my $summary = 1;
my $verbose;
my $public_type = 'FILTERED_FASTQ';
my $collection_type = 'STUDY_TYPE';
my $study_id_file;
my $update_collections;
my $log_dir;
my $ftp_root = '/nfs/1000g-archive/vol1/';
my @other_types;
my $check_status;
my $check_sample;
my $check_unidentified;

&GetOptions( 
	    'dbhost=s'      => \$dbhost,
	    'dbname=s'      => \$dbname,
	    'dbuser=s'      => \$dbuser,
	    'dbpass=s'      => \$dbpass,
	    'dbport=s'      => \$dbport,
	    'era_dbuser=s' =>\$era_dbuser,
	    'era_dbpass=s' => \$era_dbpass,
	    'run!' => \$run,
	    'all_checks!' => \$all_checks,
	    'store_new!' => \$store_new,
	    'update_existing!' => \$update_existing,
	    'summary!' => \$summary,
	    'verbose!' => \$verbose,
	    'public_type=s' => \$public_type,
	    'help!' => \$help,
	    'collection_type:s' => \$collection_type,
	    'study_id_file:s' => \$study_id_file,
	    'update_collections!' => \$update_collections,
	    'log_dir:s' => \$log_dir,
	    'other_types:s@' => \@other_types,
	    'ftp_root:s' => \$ftp_root,
	    'check_status!' => \$check_status,
	    'check_sample!' => \$check_sample,
	    'check_unidentified!' => \$check_unidentified,
	   );



if($help){
  useage();
}
my $original_run = $run;
$summary = 1 if($verbose);

if($all_checks){
  $store_new = 1;
  $update_existing = 1;
  $update_collections = 1;
  $check_sample = 1;
  $check_status = 1;
  $check_unidentified = 1;
}

if($update_existing || $store_new){
  $update_collections = 1;
}

if(($store_new + $update_existing + $update_collections + $check_sample + 
    $check_status + $check_unidentified) == 0){
  print STDERR "There script has nothing to do you need to specify at least ".
    "one of -store_new, -update_existing, -update_collections, -check_sample ".
      "-check_status or -all to have everything run\n";
  exit(0);
}




my $db = get_erapro_conn($era_dbuser, $era_dbpass);

my $reseq_db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $updater = ReseqTrack::Tools::UpdateRunMetaInfo->new
  (
   -era_db => $db,
   -dcc_db => $reseq_db,
   -study_id_file => $study_id_file,
   -log_dir => $log_dir,
   -verbose => $verbose,
   -collection_type => $collection_type,
   -filtered_fastq_type => $public_type,
   -ftp_root => $ftp_root,
   -other_types => \@other_types, 
  );


$updater->find_differences if($store_new || $update_existing);
if($store_new){
  $updater->store_new_entries;
}
if($update_existing){
  $updater->update_existing_entries;
}
if($update_collections){
  $updater->update_collections;
}
my ($dcc_problems,$era_problems);

if($check_unidentified){
  ($era_problems, $dcc_problems) = $updater->check_unidentified;
}
$updater->close_logging_fh;
my $log_file = $updater->logging_filepath;
open(FH, ">>".$log_file) or throw("Failed to open ".$log_file);
my $sample_count;
if($check_sample){
  my $move_hash =  $updater->check_sample_in_path();
  $sample_count = 0;
  $sample_count = keys(%$move_hash) if(keys(%$move_hash));
  foreach my $key(keys(%$move_hash)){
    print FH $key."\t".$move_hash->{$key}."\n";
  }
}
my $status_count;
if($check_status){
  my $move_hash = $updater->check_status();
  $status_count = 0;
  $status_count = keys(%$move_hash) if(keys(%$move_hash));
  foreach my $key(keys(%$move_hash)){
    #print $key."\t".$move_hash->{$key}."\n";
    print FH $key."\t".$move_hash->{$key}."\n";
  }
}
close(FH);
print "There are ".$sample_count." sample issues to resolve\n" if($sample_count);
print "There are ".$status_count." status issues to resolve\n" if($status_count);

print "There are ".$era_problems." sample name issues\n" if($era_problems && $verbose);
print "There are ".$dcc_problems." sample name issues\n" if($dcc_problems);

print $updater->logging_filepath."\n" if($sample_count || $status_count || 
					 $dcc_problems);

=pod

=head1 NAME

ReseqTrack/scripts/run_meta_info/update_run_meta_info_from_era.pl

=head1 SYNOPSIS

This script compares the contents of the run_meta_info table with the contents
of g1k_sequence_index from the ERAPRO database

It will summarise any new entries and differences and store them if you
let it. It can also fix sample swaps though this is mostly untested functionality.

Unless -run is specified the script checks with the user to see if a particular
update should be run

=head1 OPTIONS

-dbhost, the name of the mysql-host

-dbname, the name of the mysql database

-dbuser, the name of the mysql user

-dbpass, the database password if appropriate

-dbport, the port the mysql instance is running on, this defaults to 4197 the 
    standard port for mysql-g1kdcc.ebi.ac.uk

-era_dbuser, this the user for the ERAPRO database

-era_dbpass, this is the password for the ERAPRO database

-run, this flag means all the updates are run without first checking interactively

-store_new, this flag means new objects from the ERA database are checked and 
potentially stored

-update_existing, this flag means objects where there is a difference between the ERA
database and the DCC database are updated to match the ERA

-fix_sample, this flag means samples are checked to ensure they are correct with the
ERA database and otherwise fixed

-withdraw_suppressed, this flag means any suppressed ERA entries are checked to see
if there are associated public files, if there are they are withdrawn

-check_individual, this flag means that individuals are checked in the file paths
                   as well as at the run meta info level so filesystem differences
                   are spotted

-all_checks, this switches on all the different checks

-summary, this means summaries of the differences are printed, this is on by default

-verbose, this means details of the differences are printed, this is off by default
 switch it on also switches on the summary info

-public_type, this is the type from the filesystem which is considered when checking
              for public fastqs for suppressed runs. It defaults to FILTERED_FASTQ

-collection_type, this is the type of collection object the runs should be grouped
                  into

-study_id_list, this is the study id list which defined what makes up a collection

-help, this is a binary flag to print out the help

=head1 Examples


perl ReseqTrack/scripts/run_meta_info/update_run_meta_info_from_era.pl -dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database -era_dbuser erauser -era_dbpass ***** -store_new 


=head1 Other useful scripts

ReseqTrack/scripts/run_meta_info/dump_sequence_index.pl, This script dumps an index
file based on the current status of the file table

=cut
