#!/sw/arch/bin/perl 

use strict;
use ReseqTrack::Tools::ERAUtils;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::ERADBAdaptor;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::RunMetaInfoUtils;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::Intersection;
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
	    'fix_sample!' => \$fix_sample,
	    'withdraw_suppressed!' => \$withdraw_suppressed,
	    'check_individual!' => \$check_individual,
	    'summary!' => \$summary,
	    'verbose!' => \$verbose,
	    'public_type=s' => \$public_type,
	    'help!' => \$help,
	    'collection_type:s' => \$collection_type,
	    'study_id_file:s' => \$study_id_file,
	    'update_collections!' => \$update_collections,
    );

if($help){
  useage();
}
my $original_run = $run;
$summary = 1 if($verbose);
if($all_checks){
  $store_new = 1;
  $update_existing = 1;
  $fix_sample = 1;
  $withdraw_suppressed = 1;
  $check_individual = 1;
  $update_collections = 1;
}

if($update_existing || $store_new){
  $update_collections = 1;
}

if(($store_new + $update_existing + $fix_sample + $withdraw_suppressed + $check_individual + $update_collections) == 0){
  print STDERR "There script has nothing to do you need to specify at least ".
      "one of -store_new, -update_existing, -fix_sample, -withdraw_suppressed ".
      " or -all to have everything run\n";
  exit(0);
}

print STDERR "Trying to get a connection with ".$era_dbuser." and ".$era_dbpass."\n";
my $db = get_erapro_conn($era_dbuser, $era_dbpass);

my $reseq_db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

#Fetching meta info
my $era_rmi_a = $db->get_ERARunMetaInfoAdaptor;
print "FETCHING ALL FROM ERAPRO\n";
my $era_meta_info = $era_rmi_a->fetch_all if($store_new || $update_existing ||
					     $fix_sample);
my $era_hash = get_meta_info_hash($era_meta_info);
my $dcc_rmi_a = $reseq_db->get_RunMetaInfoAdaptor;
my $dcc_meta_info = $dcc_rmi_a->fetch_all;
my $dcc_hash = get_meta_info_hash($dcc_meta_info);

#Comparing data sets
my ($era_only, $dcc_only, $diff) = compare_era_and_dcc_meta_info($era_meta_info,
                                                                 $dcc_meta_info)
  if($store_new || $update_existing || $fix_sample);

if($summary){
  print "There are ".keys(%$era_only)." run ids only in the ERA\n";
  print "There are ".keys(%$dcc_only)." run ids only in the DCC\n";
  print "There are ".keys(%$diff)." run ids with differences between the two ".
      "databases\n";
}

my %suppressed;
foreach my $run_id(keys(%$era_hash)){
  $suppressed{$run_id} = $era_hash->{$run_id} 
  if($era_hash->{$run_id}->status eq 'suppressed');
}

#store new entries
if($store_new){
  if($summary){
    print "There are ".keys(%$era_only)." new runs to store in ".$dbname."\n";
  }
  
  foreach my $run_id (keys(%$era_only)){
    print rmi_summary_string($era_only->{$run_id})."\n" if($verbose);
  }
  
  unless($run){
      print STDERR "\n";
      print STDERR "Do you want to store these runs\n";
      print STDERR "Answer y or n :\n";
      if(get_input_arg()){
        $run = 1;
      }
    }
  if($run){
    foreach my $run_id (keys(%$era_only)){
      my $object = $era_only->{$run_id};
      $dcc_rmi_a->store($object);
    }
  }
}
$run = $original_run;

my %sample_diff;
my $print_warning;
foreach my $run_id(keys(%$diff)){
  my $era = $era_hash->{$run_id};
  my $dcc = $dcc_hash->{$run_id};
  if($era->submission_date =~ /[a-z][A-Z]/){
    throw($run_id." has a submission date ".$era->submission_date.
          " which contains letters");
  }
  if($era->sample_name ne $dcc->sample_name){
    print STDERR $era->run_id." sample name is changing from ".$dcc->sample_name.
        " \nYou must run this script with -fix_sample specific \n" unless($print_warning);
    $print_warning = 1;
    $sample_diff{$run_id} = 1;
  }
}
if($update_existing){
  if($summary){
    print "There are ".keys(%$diff)." runs with differences to account for\n";
  }
  #if($verbose){
  foreach my $run_id (keys(%$diff)){
    print rmi_summary_string($diff->{$run_id})."\n" if($verbose);
    my $dcc = $dcc_hash->{$run_id};
    my $era = $era_hash->{$run_id};
    $era->dbID($dcc->dbID);
    #if($dcc->archive_read_count >= 1 && $era->archive_read_count <= 0){
    #  throw("Trying to set ".$run_id." archive_read_count to zero");
    #}
    if($dcc->archive_base_count >= 1 && $era->archive_base_count <= 0){
      #throw("Trying to set ".$run_id." archive_base_count to zero");
      $era->archive_base_count($dcc->archive_base_count);
      $era->archive_read_count($dcc->archive_read_count);
    }
    my $history = create_history_for_run_meta_info($era, $dcc);
    if(!$history){
      throw("Can't find a difference bewteen era and dcc objects for ".$run_id);
    }else{
      $era->history($history);
      print $history->comment."\n" if($verbose);
    }
  }
  #}
  unless($run){
      print STDERR "\n";
      print STDERR "Do you want to update these runs\n";
      print STDERR "Answer y or n :\n";
      if(get_input_arg()){
        $run = 1;
      }
  }
  if($run){
    foreach my $run_id (keys(%$diff)){
      my $era = $era_hash->{$run_id};
      throw("Can't update ".$era->run_id." is has no dbID")
          unless($era->dbID);
      #print "Updating ".$era->run_id." setting base count = ".$era->archive_base_count."\n";
      $dcc_rmi_a->update($era);
    }
  }
}
$run = $original_run;

if($update_collections){
  my $ca = $reseq_db->get_CollectionAdaptor;
  my $study_id_hash = parse_study_id_list($study_id_file);
  my $rmis = $dcc_rmi_a->fetch_all;
  my %study_id_to_rmi;
  foreach my $rmi(@$rmis){
    $study_id_to_rmi{$rmi->study_id}{$rmi->run_id} = $rmi;
  }
  my $existing_collections = $ca->fetch_by_type($collection_type);
  my %collection_hash;
  foreach my $collection(@$existing_collections){
    $collection_hash{$collection->name} = $collection;
  }
  unless($run){
    print STDERR "\n";
    print STDERR "Do you want to update the collections on these runs\n";
    print STDERR "Answer y or n :\n";
    if(get_input_arg()){
      $run = 1;
    }else{
      exit(0);
    }
  }
  foreach my $name(keys(%$study_id_hash)){
    my $collection = $collection_hash{$name};
    my $study_list_hash = $study_id_hash->{$name};
    my @rmis;
    my %run_id_to_rmi;
    foreach my $study_id(keys(%$study_list_hash)){
      my $rmi_hash = $study_id_to_rmi{$study_id};
      my @rmis_from_study = values(%$rmi_hash);
      foreach my $rmi(@rmis_from_study){
	$run_id_to_rmi{$rmi->run_id} = $rmi;
      }
      push(@rmis, @rmis_from_study);
    }
    unless($collection){
      $collection = ReseqTrack::Collection->new(
						-name => $name,
						-others => \@rmis,
						-type => $collection_type,
						-table_name => 'run_meta_info',
					       );
      $collection_hash{$name} = $collection;
    }else{
      my $collection_rmis = $collection->others;
      my @collection_run_ids;
      my @run_ids = keys(%run_id_to_rmi);
      foreach my $other(@$collection_rmis){
	push(@collection_run_ids, $other->run_id);
      }
      my $file_set = ReseqTrack::Tools::Intersection->new(
							  -list => \@run_ids,
							 );
      my $db_set = ReseqTrack::Tools::Intersection->new
	(
	 -list => \@collection_run_ids,
	);
      my $not_in_file = $db_set->not($file_set);
      my $not_in_db = $file_set->not($db_set);
      my $db_and_file = $db_set->and($file_set);
      #add new to collection
      foreach my $run_id(@{$not_in_db->list}){
	my $rmi = $run_id_to_rmi{$run_id};
	throw("Failed to find a rmi for ".$run_id) unless($rmi);
	$collection->other($collection);
      }
      $collection_hash{$name} = $collection;
      #remove old from collection
      my @rmis_to_remove;
      foreach my $run_id(@{$not_in_file->list}){
	my $rmi = $run_id_to_rmi{$run_id};
	push(@rmis_to_remove, $rmi);
      }
      $ca->remove_others($collection, \@rmis_to_remove) if($run);
    }
  }
  foreach my $name(keys(%collection_hash)){
    my $collection = $collection_hash{$name};
    $ca->store($collection, 1);
  }
}
$run = $original_run;
#fix sample swaps
if($fix_sample){
  #print summary
  #ask if want to fix
  print STDERR "*****fix sample THIS IS TOTALLY UNTESTED CODE be certain that you want to run it ******\n";
  if($summary){
    print "There are ".keys(%sample_diff)." runs with sample differences which ".
        "need to be fixed\n";
    if($verbose){
      my @keys = keys(%sample_diff);
      foreach my $run_id(keys(%sample_diff)){
        my $era = $era_hash->{$run_id};
        my $dcc = $dcc_hash->{$run_id};
        print $run_id." era sample ".$era->sample_name." dcc sample ".$dcc->sample_name." ".$era->center_name." ".$era->study_id."\n";
      }
    }
  }
  unless($run){
      print STDERR "\n";
      print STDERR "Do you want to fix the sample on these runs\n";
      print STDERR "Answer y or n :\n";
      if(get_input_arg()){
        $run = 1;
      }else{
	exit(0);
      }
  }
  if($run){
    foreach my $run_id(keys(%sample_diff)){
      my $era = $era_hash->{$run_id};
      my $dcc = $dcc_hash->{$run_id};
      fix_sample_swap($run_id, $dcc->sample_name, $era->sample_name, $db, $run);
    }
  }
}
$run = $original_run;
#withdraw suppressed files

$dcc_meta_info = $dcc_rmi_a->fetch_all;
my $fa =  $reseq_db->get_FileAdaptor;
my $ca = $reseq_db->get_CollectionAdaptor;
if($withdraw_suppressed){
  print STDERR "***** withdraw suppressed THIS IS TOTALLY UNTESTED CODE be certain that you want to run it ******\n";
  my @to_withdraw;
  foreach my $meta_info(@$dcc_meta_info){
    my $run_id  = $meta_info->run_id;
    if($meta_info->status eq 'suppressed'){
      my $files = $fa->fetch_all_like_name($run_id);
      if($files && @$files >= 1){
        foreach my$file(@$files){
          next unless($file->type eq $public_type);
          next if($file->name =~ /withdrawn/);
          push(@to_withdraw, $file);
        }
      }
    }
  }
  if(@to_withdraw >= 1){
    print "There are ".@to_withdraw." files which need to be withdrawn\n" if($summary);
    if($verbose){
      foreach my $file(@to_withdraw){
        print $file->name." needs to be withdrawn\n";
      }
      print "Need to implement automated withdrawal\n";
    }
  }
}
$run = $original_run;
if($check_individual){
  print STDERR "*****check individual THIS IS TOTALLY UNTESTED CODE be certain that you want to run it ******\n";
  my %name_to_individual;
  my @files_to_fix;
  my $archive_fastq = $fa->fetch_by_type("ARCHIVE_FASTQ");
  my $filtered_fastq = $fa->fetch_by_type("FILTERED_FASTQ");
  my @files;
  push(@files, @$archive_fastq);
  push(@files, @$filtered_fastq);
  my %run_to_files;
  foreach my $file(@files){
    my $filename = $file->filename;
    $filename =~ /^([E|S]RR\d+)\./;
    my $run_id = $1;
    push(@{$run_to_files{$run_id}}, $file);
  }
  foreach my $meta_info(@$dcc_meta_info){
    next if($meta_info->status eq 'suppressed');
    my $run_id = $meta_info->run_id;
    my $sample_name = $meta_info->sample_name;
    my $files =  $run_to_files{$run_id};
    if($files && @$files >= 1){
      foreach my $file(@$files){
        next unless($file->type eq $public_type);
        next if($file->name =~ /withdrawn/);
	if($file->name =~ /$sample_name/){
	  next;
	}
        $name_to_individual{$file->name} = $sample_name;
        push(@files_to_fix, $file);
      }
    }
  }
  print "There are ".@files_to_fix." files with the wrong individual in their ".
      "path\n" if($summary);
  if($verbose){
    foreach my $file(@files_to_fix){
      print $file->name." should be in ".$name_to_individual{$file->name}."'s dir\n";
    }
  }
  unless($run){
      print STDERR "\n";
      print STDERR "Do you want to fix the sample directories on these runs\n";
      print STDERR "Answer y or n :\n";
      if(get_input_arg()){
        $run = 1;
      }
  }
  if($run){
    my %run_ids;
    foreach my $file(@files_to_fix){
      $file->name =~ /([E|S]RR\d+)/;
      my $run_id = $1;
      next if($run_ids{$run_id});
      $run_ids{$run_id} = 1;
      $file->name =~ /(NA\d+)/;
      my $old_sample = $1;
      my $new_sample = $name_to_individual{$file->name};
      #fix_sample_swap($run_id, $old_sample, $new_sample, $db, $run);
    }
  }
}

if(keys(%{$dcc_only}) >= 1){
  print STDERR "There is a problem there are runs which are only in the dcc ".
      "database ".$dbname."\n";
  exit(1);
}

sub rmi_summary_string{
  my ($object) = @_;
  my $string =  $object->run_id." ".$object->submission_id." ".$object->center_name." ".$object->study_id;
  return $string;
}
sub useage{
  exec('perldoc', $0);
  exit(0);
}

sub parse_study_id_list{
  my ($study_id_list) = @_;
  open(FH, "<", $study_id_list) or throw("Failed to open ".$study_id_list." $!");
  my %config;
  my %seen;
  my $header;
  while(<FH>){
    chomp;
    next if (/^\s$/ || /^\#/);
    if (/^\[(.*)\]\s*$/) {    # $1 will be the header name, without the []
      $header = $1;
      if ($config{$header}) {
	throw("Shouldn't have duplicate headers in ".$study_id_list." ".$header);
      }
      $config{$header} = {};
    }else{
      my $line = trim_spaces($_);
      if($line){
	if($seen{$line}){
	  throw($line." seems to be duplicated in ".$study_id_list);
	}else{
	  $seen{$line} = 1;
	}
	$config{$header}{$line} = 1;
      }else{
	warning($header." has blank lines");
      }
    }
  }
  return \%config;
}

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
