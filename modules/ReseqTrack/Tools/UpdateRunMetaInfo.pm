=head1 NAME

ReseqTrack::Tools::UpdateRunMetaInfo

=head1 SYNOPSIS

Object to update and sanity check the run meta info tables

=head1 Example


=cut

package ReseqTrack::Tools::UpdateRunMetaInfo;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::ERAUtils;
use ReseqTrack::Tools::RunMetaInfoUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::FileUtils;
use File::Path;
use FileHandle;
use Data::Dumper;

sub new{
  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;
  my ($era_db, $dcc_db, $era_dbuser, $era_dbpass, $dbname, $dbuser, $dbport, $dbpass,
      $dbhost, $store_new, $update, $study_id_file, $log_dir, $verbose, 
      $collection_type, $filtered_fastq_type, $ftp_root, $other_types,
      $population_rules, $study_ids) = 
	rearrange([qw(ERA_DB DCC_DB ERA_DBUSER ERA_DBPASS DBNAME DBUSER DBPORT 
		      DBPASS DBHOST STORE_NEW UPDATE STUDY_ID_FILE LOG_DIR VERBOSE 
		      COLLECTION_TYPE FILTERED_FASTQ_TYPE FTP_ROOT OTHER_TYPES
                      POPULATION_RULES STUDY_IDS)], 
		  @args);
  #SETTING DEFAULTS
  $self->filtered_fastq_type('FILTERED_FASTQ');
  $self->other_types(['ARCHIVE_FASTQ', 'WITHDRAWN_FILTERED_FASTQ']);
  $self->collection_type('STUDY_TYPE');
  $self->ftp_root('/nfs/1000g-archive/vol1/');
  #######
  $self->era_dbuser($era_dbuser);
  $self->era_dbpass($era_dbpass);
  $self->era_db($era_db);
  $self->dbhost($dbhost);
  $self->dbname($dbname);
  $self->dbuser($dbuser);
  $self->dbpass($dbpass);
  $self->dbport($dbport);
  $self->dcc_db($dcc_db);
  $self->verbose($verbose);
  $self->study_id_file($study_id_file);
  $self->log_dir($log_dir);
  $self->collection_type($collection_type);
  $self->filtered_fastq_type($filtered_fastq_type);
  $self->ftp_root($ftp_root);
  $self->population_rules($population_rules) if $population_rules;
  $self->study_ids($study_ids) if $study_ids;
  
  return $self;
}

sub study_id_file{
  my ($self, $arg) = @_;
  if($arg){
    $self->{study_id_file} = $arg;
  }
  return $self->{study_id_file};
}

sub dcc_db {
 my ( $self, $arg ) = @_;
 if($arg){
   $self->{dcc_db} = $arg;
 }
 unless($self->{dcc_db}){
   $self->{dcc_db} = ReseqTrack::DBSQL::DBAdaptor->new(
						       -host   => $self->dbhost,
						       -user   => $self->dbuser,
						       -port   => $self->dbport,
						       -dbname => $self->dbname,
						       -pass   => $self->dbpass,
						      );
 }
}


sub ftp_root{
  my ($self, $arg) = @_;
  $self->{ftp_root} = $arg if($arg);
  return $self->{ftp_root};
}

sub dbhost {
 my ( $self, $arg ) = @_;
   $self->{dbhost} = $arg if ( $arg);
 return $self->{dbhost};
}

sub dbname {
 my ( $self, $arg ) = @_;
   $self->{dbname} = $arg if ( $arg);
  return $self->{dbname};
}

sub dbuser {
 my ( $self, $arg ) = @_;
   $self->{dbuser} = $arg if ( $arg);
  return $self->{dbuser};
}

sub dbpass {
 my ( $self, $arg ) = @_;
  $self->{dbpass} = $arg if ( $arg);
 return $self->{dbpass};
}

sub dbport {
 my ( $self, $arg ) = @_;
  $self->{dbport} = $arg if ($arg);
 return $self->{dbport};
}

sub era_dbuser{
  my ( $self, $arg ) = @_;
  $self->{era_dbuser} = $arg if ($arg);
  return $self->{era_dbuser};
}

sub era_dbpass{
  my ( $self, $arg ) = @_;
  $self->{era_dbpass} = $arg if ($arg);
  return $self->{era_dbpass};
}

sub era_db{
  my ($self, $arg) = @_;
  if($arg){
    $self->{era_db} = $arg;
  }
  unless($self->{era_db}){
    $self->{era_db} =  get_erapro_conn($self->era_dbuser, $self->era_dbpass);
  }
}

sub log_dir{
  my ($self, $arg) = @_;
  if($arg){
    mkpath($arg) unless(-d $arg);
    $self->{log_dir} = $arg;
    throw($arg." must be a directory") unless(-d $arg);
  }
  return $self->{log_dir};
}


sub verbose{
  my ($self, $arg) = @_;
  if(defined($arg)){
    $self->{verbose} = $arg;
  }
  return $self->{verbose};
}

sub logging_fh{
  my ($self, $arg) = @_;
  if($arg){
    $self->{logging_fh} = $arg;
  }
  unless($self->{logging_fh}){
    my $open_string = ">".$self->logging_filepath;
    my $fh = FileHandle->new($open_string);
    my ($p, $f, $l) = caller;
    unless($fh){
      throw("Failed to open ".$self->logging_filepath." $!");
    }
    $self->{logging_fh} = $fh;
  }
  return $self->{logging_fh};
}

sub close_logging_fh{
  my ($self) = @_;
  $self->{logging_fh} = undef;
}

sub logging_filepath{
  my ($self, $arg) = @_;
  if($arg){
    $self->{logging_filepath} = $arg;
  }
  unless($self->{logging_filepath}){
    my $ident = int(rand(10000));
    my $date = current_date();
    my $tmp_name = "update_runmetainfo.".$date.".".$ident.".$$.log";
    $self->{logging_filepath} = $self->log_dir."/".$tmp_name;
  }
  return $self->{logging_filepath};
}

sub era_only{
  my ($self, $arg) = @_;
  if($arg){
    $self->{era_only} = $arg;
  }
  return $self->{era_only};
}

sub dcc_only{
  my ($self, $arg) = @_;
  if($arg){
    $self->{dcc_only} = $arg;
  }
  return $self->{dcc_only};
}

sub diff_set{
  my ($self, $arg) = @_;
  if($arg){
    $self->{diff_set} = $arg;
  }
  return $self->{diff_set};
}

sub population_rules {
  my ($self, $rules, $force_update) = @_;
  if ($rules) {
    throw("rules must be an arrayref") if (ref $rules ne 'ARRAY');
    $self->{population_rules} = $rules;
  }
  if (!$self->{rules} || $force_update){
    $self->{population_rules} = $self->dcc_db->get_PopulationRuleAdaptor->fetch_all_in_order();
  }
  return $self->{population_rules};
}

sub study_ids {
  my ($self, $study_ids, $force_update) = @_;
  if ($study_ids) {
    throw("study_ids must be an arrayref") if (ref $study_ids ne 'ARRAY');
    $self->{study_ids} = $study_ids;
  }
  if(!$self->{study_ids} || $force_update){
    $self->{study_ids} = $self->dcc_db->get_StudyIDAdaptor->fetch_all();
  }
  return $self->{study_ids};
}

sub find_differences{
  my ($self) = @_;
  my $era_rmis = $self->get_era_rmis();
  my $dcc_rmis = $self->get_dcc_rmis();
  my $logging_fh = $self->logging_fh();
  my ($era_only, $dcc_only, $diff) = compare_era_and_dcc_meta_info($era_rmis,
								   $dcc_rmis);
  #if($self->verbose){
  print $logging_fh "There are ".keys(%$era_only)." runs only in the ERA ".
    "database\n";
  print $logging_fh "There are ".keys(%$dcc_only)." runs only in the DCC ".
    "database\n";
  print $logging_fh "There are ".keys(%$diff)." differences in runs in both ".
    "the ERA and the DCC database\n";
  #}
  $self->era_only($era_only);
  $self->dcc_only($dcc_only);
  $self->diff_set($diff);
}

sub sample_diffs{
  my ($self, $arg) = @_;
  if($arg){
    $self->{sample_diffs} = $arg;
  }
  unless($self->{sample_diffs}){
    $self->{sample_diffs} = $self->find_sample_diffs;
  }
  return $self->{sample_diffs};
}

sub store_new_entries{
  my ($self, $rmis) = @_;
  $rmis = $self->era_only unless($rmis);
  my $rmia = $self->dcc_db->get_RunMetaInfoAdaptor;
  my $logging_fh = $self->logging_fh();
  foreach my $run_id(keys(%$rmis)){
    my $rmi = $rmis->{$run_id};
    print $logging_fh "STORING ".$self->rmi_summary_string($rmi)."\n" if($self->verbose);
    $rmia->store($rmi);
  }
  return 1;
}

sub update_existing_entries{
  my ($self, $rmis) = @_;
  $rmis = $self->diff_set unless($rmis);
  my $dcc_hash = $self->dcc_hash;
  my $rmia = $self->dcc_db->get_RunMetaInfoAdaptor;
  my $logging_fh = $self->logging_fh();
  foreach my $run_id(keys(%$rmis)){
    my $rmi = $rmis->{$run_id};
    print $logging_fh "DIFFS in ".$self->rmi_summary_string($rmi)."\n" if($self->verbose);
    my $dcc = $dcc_hash->{$rmi->run_id};
    if($dcc->archive_base_count >= 1 && $rmi->archive_base_count <= 0){
      $rmi->archive_base_count($dcc->archive_base_count);
    }
    if($rmi->sample_name eq 'unidentified'){
      $rmi->sample_name($dcc->sample_name);
    }
    $rmi->dbID($dcc->dbID);
    my $history = create_history_for_run_meta_info($rmi, $dcc);
    next unless($history);
    $rmi->history($history);
    print $logging_fh $rmi->run_id." has changed ".$history->comment."\n" if($self->verbose);
    $rmia->update($rmi);
  }
}

sub update_collections{
  my ($self, $study_id_file) = @_;
  $study_id_file = $self->study_id_file unless($study_id_file);
  my $study_id_hash = $self->parse_study_id_file($study_id_file);
  my $ca = $self->dcc_db->get_CollectionAdaptor;
  my $rmis = $self->get_dcc_rmis(undef, 1);
  my %study_id_to_rmi;
  foreach my $rmi(@$rmis){
    $study_id_to_rmi{$rmi->study_id}{$rmi->run_id} = $rmi;
  }
  my $existing_collections = $ca->fetch_by_type($self->collection_type);
  my %collection_hash;
  foreach my $collection(@$existing_collections){
    $collection_hash{$collection->name} = $collection;
  }
  foreach my $name(keys(%$study_id_hash)){
    my $collection = $collection_hash{$name};
    my $study_id_list = $study_id_hash->{$name};
    my @rmis;
    my %run_id_to_rmi;
    foreach my $study_id(keys(%$study_id_list)){
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
						-type => $self->collection_type,
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
	$collection->others($rmi);
      }
      $collection_hash{$name} = $collection;
      #remove old from collection
      my @rmis_to_remove;
      foreach my $run_id(@{$not_in_file->list}){
	my $rmi = $run_id_to_rmi{$run_id};
	push(@rmis_to_remove, $rmi) if($rmi);
      }
      $ca->remove_others($collection, \@rmis_to_remove) if(@rmis_to_remove >= 1);
    }
  }
  foreach my $name(keys(%collection_hash)){
    my $collection = $collection_hash{$name};
    $ca->store($collection, 1);
  }
}

sub check_status{
  my ($self, $rmis, $files) = @_;
  $rmis = $self->get_dcc_rmis() unless($rmis);
  $files = $self->get_filtered_fastq unless($files);
  my $file_hash = $self->run_hash_from_fastq($files) if($files && @$files);
  my $ftp_root = $self->ftp_root;
  my %move_hash;
  foreach my $rmi(@$rmis){
    next if($rmi->status eq 'public');
    my $files = $file_hash->{$rmi->run_id};
    foreach my $file(@$files){
      next unless($file->name =~ /$ftp_root/);
      next if($file->name =~ /withdrawn/);
      #print "Looking at ".$file->name." ".$file->dbID." ".$file->type."\n";
      my $new = $file->name;
      $new =~ s/ftp/withdrawn/;
      next if($new eq $file->name);
      $move_hash{$file->name} = $new;
    }
  }
  return \%move_hash;
}

sub check_sample_in_path{
  my ($self, $rmis, $files) = @_;
  $rmis = $self->get_dcc_rmis() unless($rmis);
  $files = $self->get_all_fastq unless($files) ;
  return {} unless($files && @$files >= 1);
  my $file_hash = $self->run_hash_from_fastq($files) if($files && @$files);
  my %move_hash;
  my $verbose = 0;
  foreach my $rmi(@$rmis){
    my @problem_files;
    my $sample_name = $rmi->sample_name;
    my $files = $file_hash->{$rmi->run_id};
    foreach my $file(@$files){
      next if($file->name =~ /$sample_name/);
      #print "Fail to find ".$sample_name." in ".$file->name."\n";
      push(@problem_files, $file);
    }
    foreach my $file(@problem_files){
      my $name = $file->name;
      my $new_name = $file->name;
      $name =~ /(NA\d+)/;
      my $old_sample = $1;
      unless($old_sample){
	$name =~ /(HG\d+)/;
	$old_sample = $1;
      }
      unless($old_sample){
	if($name =~ /unidentified/){
	  $old_sample = 'unidentified';
	}
      }
      if(!$old_sample){
	throw("Failed to get a sample id from ".$file->name);
      }
      $new_name =~ s/$old_sample/$sample_name/;
      $new_name =~ s/withdrawn/ftp/;
      if($new_name =~ /withdrawn/){
	throw("Still have withdrawn in new name ".$new_name);
      }
      if($name eq $new_name){
	throw("There is a problem ".$name." is the same as ".$new_name." when ".
	      "changing from ".$old_sample." to ".$sample_name);
      }
      $move_hash{$name} = $new_name;
    }
  }
  return \%move_hash;
}

sub sort_move_hash{
  my ($self, $hash) = @_;
  my $ftp_root = $self->ftp_root;
  my %ftp;
  my %non_ftp;
  foreach my $name(keys(%$hash)){
    if($name =~ /$ftp_root/){
      $ftp{$name} = $hash->{$name};
    }else{
      $non_ftp{$name} = $hash->{$name};
    }
  }
  return (\%ftp, \%non_ftp);
}

sub fix_file_type{
  my ($self, $files) = @_;
  $files = $self->get_all_fastq unless($files);
  my $fa = $self->dcc_db->get_FileAdaptor;
  foreach my $file(@$files){
    if($file->type =~ /withdrawn/i){
      next if($file->name =~ /withdrawn/);
      my $new_file = copy_file_object($file);
      my $new_type = $file->type;
      $new_type =~ s/WITHDRAWN\_//g;
      $new_file->type($new_type);
      my $history = create_history($new_file,$file);
      next unless($history);
      $new_file->history($history);
      $fa->update($new_file);
    }
  }
}

sub print_move_hash{
  my ($self, $hash, $fh) = @_;
  my ($p, $f, $l) = caller;
  print $f.":".$l." calling print move hash\n";
  $fh = \*STDOUT unless($fh);
  foreach my $key(keys(%$hash)){
    print $fh $key."\t".$hash->{$key}."\n";
  }
}

sub get_filtered_fastq{
  my ($self, $files) = @_;
  if($files){
    $self->{filtered_fastq} = $files
  }
  unless($self->{filtered_fastq}){
    my $fa = $self->dcc_db->get_FileAdaptor;
    $self->{filtered_fastq} = $fa->fetch_by_type($self->filtered_fastq_type);
  }
}

sub get_all_fastq{
  my ($self, $files) = @_;
  if($files){
    $self->{all_fastq} = $files;
  }
  unless($self->{all_fastq}){
    my @files;
    push(@files, @{$self->get_filtered_fastq});
    my $fa = $self->dcc_db->get_FileAdaptor;
    foreach my $type(@{$self->other_types}){
      push(@files, @{$fa->fetch_by_type($type)});
    }
    $self->{all_fastq} = \@files;
    #throw("Have no files in database ") unless(@{$self->{all_fastq}});
  }
  return $self->{all_fastq};
}

sub filtered_fastq_type{
  my ($self, $type) = @_;
  if($type){
    $self->{filtered_fastq_type} = $type;
  }
  return $self->{filtered_fastq_type};
}

sub other_types{
  my ($self, $type) = @_;
  if($type){
    $self->{other_types} = $type;
  }
  return $self->{other_types};
}

sub collection_type{
  my ($self, $type) = @_;
  if($type){
    $self->{collection_type} = $type;
  }
  return $self->{collection_type};
}
sub parse_study_id_file{
  my ($self, $file, $reparse) = @_;
  if(!$self->{study_id_hash} || $reparse){
    $file = $self->study_id_file unless($file);
    open(FH, "<", $file) or throw("Failed to open ".$file." $!");
    my %config;
    my %seen;
    my $header;
    while(<FH>){
      chomp;
      next if (/^\s$/ || /^\#/);
      if (/^\[(.*)\]\s*$/) {    # $1 will be the header name, without the []
	$header = $1;
	if ($config{$header}) {
	  throw("Shouldn't have duplicate headers in ".$file." ".$header);
	}
	$config{$header} = {};
      }else{
	my $line = trim_spaces($_);
	if($line){
	  if($seen{$line}){
	    throw($line." seems to be duplicated in ".$file);
	  }else{
	    $seen{$line} = 1;
	  }
	  $config{$header}{$line} = 1;
	}else{
	  #warning($header." has blank lines");
	}
      }
    }
    $self->{study_id_hash} = \%config;
  }
  return $self->{study_id_hash};
}

sub find_sample_diffs{
  my ($self) = @_;
  my %sample_diffs;
  my $logging_fh = $self->logging_fh();
  foreach my $run_id(keys(%{$self->diff_set})){
    my $era_rmi = $self->era_hash->{$run_id};
    my $dcc_rmi = $self->dcc_hash->{$run_id};
    unless($era_rmi->sample_name eq $dcc_rmi->sample_name){
      $sample_diffs{$run_id} = 1;
      print $logging_fh $run_id." has a sample difference ".
	$dcc_rmi->sample_name." has become ".$era_rmi->sample_name." in the ".
	  "archive\n";
    }
  }
  return \%sample_diffs;
}

sub era_hash{
  my ($self, $arg) = @_;
  if($arg){
    $self->{era_hash} = $arg;
  }
  unless($self->{era_hash}){
    my %hash;
    foreach my $rmi(@{$self->get_era_rmis()}){
      $hash{$rmi->run_id} = $rmi;
    }
    $self->{era_hash} = \%hash;
  }
  return $self->{era_hash};
}

sub dcc_hash{
  my ($self, $arg) = @_;
  if($arg){
    $self->{dcc_hash} = $arg;
  }
  unless($self->{dcc_hash}){
    my %hash;
    foreach my $rmi(@{$self->get_dcc_rmis()}){
      $hash{$rmi->run_id} = $rmi;
    }
    $self->{dcc_hash} = \%hash;
  }
  return $self->{dcc_hash};
}

sub suppressed_hash{
  my ($self, $hash) = @_;
  if($hash){
    $self->{suppressed_hash} = $hash;
  }
  unless($self->{suppressed_hash}){
    my $rmis = $self->get_era_rmis();
    my %hash;
    foreach my $rmi(@$rmis){
      $hash{$rmi->run_id} = $rmi->run_id if($rmi->status ne 'public');
    }
    $self->{suppressed_hash} = \%hash;
  }
  return $self->{suppressed_hash};
}

sub get_era_rmis{
  my ($self, $rmis, $force_update) = @_;
  if($rmis){
    $self->{era_rmis} = $rmis;
  }
  if(!$self->{era_rmis} || $force_update){
    my $era_rmia = $self->era_db->get_ERARunMetaInfoAdaptor(
                    -study_ids => $self->study_ids,
                    -population_rules => $self->population_rules);
    my $rmis = $era_rmia->fetch_all();
    my @array;
    foreach my $rmi(@$rmis){
      next unless($rmi);
      next if($rmi->sample_name =~ /genomic/i);
      next if($rmi->sample_name =~ /\./);
      push(@array, $rmi);
    }
    $self->{era_rmis} = \@array;
  }
  return $self->{era_rmis};
}

sub get_dcc_rmis{
  my ($self, $rmis, $force_update) = @_;
  if($rmis){
    $self->{dcc_rmis} = $rmis;
  }
  if(!$self->{dcc_rmis} || $force_update){
    $self->{dcc_rmis} = $self->dcc_db->get_RunMetaInfoAdaptor->fetch_all();
  }
  return $self->{dcc_rmis};
}


sub rmi_summary_string{
  my ($self, $object) = @_;
  my $string =  $object->run_id." ".$object->submission_id." ".$object->center_name." ".$object->study_id;
  return $string;
}


sub run_hash_from_fastq{
  my ($self, $files) = @_;
  my %hash;
  my ($p, $f, $l) = caller;
  foreach my $file(@$files){
    $file->filename =~ /([E|S]RR\d+)/;
    push(@{$hash{$1}}, $file);
  }
  return \%hash;
}


sub check_unidentified{
  my ($self) = @_;
  my $era_rmis = $self->get_era_rmis(undef, 1);
  my $dcc_rmis = $self->get_dcc_rmis(undef, 1);
  my @era_problems;
  my @dcc_problems;

  foreach my $rmi(@$era_rmis){
    next if($rmi->sample_name =~ /^NA/ || $rmi->sample_name =~ /^HG/);
    push(@era_problems, $rmi);
  }
  foreach my $rmi(@$dcc_rmis){
    next if($rmi->sample_name =~ /^NA/ || $rmi->sample_name =~ /^HG/);
    push(@dcc_problems, $rmi);
  }
  my $logging_fh = $self->logging_fh;
  print $logging_fh "There are ".@era_problems." era runs with unidentified samples\n" if(@era_problems);
  print "There are ".@dcc_problems." dcc runs with unidentified samples\n" if(@dcc_problems);
  print $logging_fh "There are ".@dcc_problems." dcc runs with unidentified samples\n" if(@dcc_problems);

 
  foreach my $rmi(@era_problems){
    print $logging_fh $rmi->run_id." ".$rmi->center_name." ".$rmi->study_name." ".$rmi->sample_name." ".$rmi->population."\n";
  }
 
  foreach my $rmi(@dcc_problems){
    print $logging_fh $rmi->run_id." ".$rmi->center_name." ".$rmi->study_name." ".$rmi->sample_name." ".$rmi->population."\n";
  }
  
  my $dcc_problem_count = scalar(@dcc_problems);
  my $era_problem_count = scalar(@era_problems);


  return $era_problem_count, $dcc_problem_count;
}

1;
