=pod

=head1 NAME

ReseqTrack::Tools::FileUtils;

=head1 SYNOPSIS

A method collection class which contains useful utility methods for ReseqTrack::File
object handling.

You can either import all the contained method into your name space

use ReseqTrack::Tools::FileUtils;

or just specific methods

use ReseqTrack::Tools::FileUtils qw(are_files_identical)

=cut

package ReseqTrack::Tools::FileUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::File;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::History;
use File::Copy;
use File::Basename;
use File::Path;
use File::Find ();
use ReseqTrack::Tools::GeneralUtils;
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(are_files_identical 
             create_objects_from_path_list 
             create_object_from_path
             print_file
             calculate_comment
             create_history
             withdrawn_file
             copy_file_object
             assign_type
             check_type
	     	 move_file_in_db_and_dir
             get_count_stats
	     	 get_run_id_to_file_hash
	     	 write_log);

=head2 are_files_identical

  Arg [1]   : ReseqTrack::File
  Arg [2]   : ReseqTrack::File
  Function  : compare two objects return 0 if they aren't the same and 1 if they are
  the only aspects which aren't compared are updated and history
  Returntype: 0/1
  Exceptions: none
  Example   : if(are_files_identical($file, $other_file){
     print $file." is the same as ".$other_file."\n";
   }

=cut

sub are_files_identical{
  my ($one, $two) = @_;
  throw("Must pass are_files_identical two file objects") unless($one && $two);
  #throw("Must have the same filename".$one->filename." ".$two->filename) unless($one->filename eq $two->filename);
  return 0  unless($one->filename eq $two->filename);
  return 0 unless($one->dirname eq $two->dirname);
  return 0 unless($one->path eq $two->path);
  return 0 unless($one->full_path eq $two->full_path);
  if($one->md5 && $two->md5){
      unless($one->md5 eq $two->md5){
        return 0;
    }
  }
  if($two->md5 && !($one->md5)){
    return 0;
  }
  if($one->md5 && !($two->md5)){
    warning("You are going to undefined the md5 of ".$two->name);
    return 0;
  }
  if( ($one->size || $one->size == 0) && ($two->size || $two->size ==0)){
    return 0 unless($one->size eq $two->size);
  }
  #  print "obj1 " . $one->name . " size " . $one->size . " obj2 " .  $two->name . " size " . $two->size ."\n";  
  return 0 unless($one->type eq $two->type);
  return 0 unless($one->host->name eq $two->host->name);
  return 0 unless($one->withdrawn eq $two->withdrawn);
  return 1;
}


=head2 create_objects_from_path_list

  Arg [1]   : Arrayref of directory paths
  Arg [2]   : string, file type
  Arg [3]   : ReseqTrack::Host 
  Function  : creates a ReseqTrack::File object for each path in the list
  Returntype: arrayref of ReseqTrack::File objects
  Exceptions: throws if not passed a type or a host
  Example   : my $files = create_objects_from_path_list($paths, "FASTQ", $host);

=cut


sub create_objects_from_path_list{
  my ($paths, $type, $host) = @_;
  throw("Must pass create_objects_from_path_list a ReseqTrack::Host object ".
        $host) if(!$host || !($host->isa("ReseqTrack::Host")));
  throw("Must pass create_objects_from_path_list a type string")
      unless($type);
  my @files;
  foreach my $path(@$paths){
    my $size = (-s $path);
    $size = 0 unless($size);
    my $file = ReseqTrack::File->new(
      -name => $path,
      -type => $type,
      -size => $size,
      -host => $host,
        );
    push(@files, $file);
  }
  return \@files,
}

=head2 create_object_from_path

  Arg [1]   : string, filepath
  Arg [2]   : string, file type
  Arg [3]   : ReseqTrack::Host 
  Function  : creates a ReseqTrack::File object for each path in the list
  Returntype: arrayref of ReseqTrack::File objects
  Exceptions: throws if not passed a type or a host
  Example   : my $files = create_objects_from_path_list($paths, "FASTQ", $host);

=cut


sub create_object_from_path{
  my ($path, $type, $host) = @_;
  throw("Must pass create_objects_from_path_list a ReseqTrack::Host object ".
        $host) if(!$host || !($host->isa("ReseqTrack::Host")));
  throw("Must pass create_objects_from_path_list a type string")
      unless($type);
  my $size = (-s $path);
  $size = 0 unless($size);
  my $file = ReseqTrack::File->new(
    -name => $path,
    -type => $type,
    -size => $size,
    -host => $host,
      );
  return $file;
}

=head2 calculate_comment

  Arg [1]   : ReseqTrack::File
  Arg [2]   : ReseqTrack::File
  Function  : checks for differences between two objects and returns string
  describing difference
  Returntype: string/undef
  Exceptions: none
  Example   : my $comment = calculate_comment($new_file, $old_file);

=cut



sub calculate_comment{
  my ($new, $old) = @_;
  #throw("New and Old files need to have the same name for this to be relavant".
  #      " not new ".$new->name." and old ".$old->name) if($new->name ne $old->name);
  my $result = are_files_identical($new, $old);
  unless($result == 1){
    if($new->full_path ne $old->full_path){
      return "File has changed location from ".$old->full_path." to ".$new->full_path;
    }  
    elsif(($new->md5 && $old->md5) && ($new->md5 ne $old->md5)){
      return "File md5 has changed from ".$old->md5." to ".$new->md5;
    }  
    elsif($new->size != $old->size){
      return "File has changed size from ".$new->size." to ".$old->size;
    }
    elsif($new->host->name ne $new->host->name){
      return "File has changed hosts from ".$old->host->name." to ".$new->host->name;
    }
    elsif($new->type ne $old->type){
      return "File has changed types from ".$old->type." to ".$new->type;
    }
  }
  if($old->md5 && !$new->md5){
    return "File is having md5 undefined from ".$old->md5;
  }
  if(!$old->md5 && $new->md5){
    return "File is having md5 set to ".$new->md5;
  }
  return undef;
}



=head2 create_history

  Arg [1]   : ReseqTrack::File
  Arg [2]   : ReseqTrack::File
  Function  : compares two given objects, creates a comment based on the difference
  and a history object from the comment
  Returntype: ReseqTrack::History
  Exceptions: n/a
  Example   : my $history = create_history($new, $old);
  $new->history($history);

=cut



sub create_history{
  my ($new, $old, $comment) = @_;
  throw("Cant create a history object if old ".$old." does not have a dbID") unless($old->dbID); 
  $comment = calculate_comment($new, $old) unless($comment);
  if($comment){
    my $history = ReseqTrack::History->new(
      -other_id => $old->dbID,
      -table_name => 'file',
      -comment => $comment,
        );
    return $history;
  }else{
    print STDERR "Have no comment, can't find a diff bewteen new and old\n";
  }
  return undef;
}


=head2 print_file

  Arg [1]   : ReseqTrack::File
  Function  : prints a string about the file containing, name, size, md5 and type
  Returntype: n/a
  Exceptions: n/a
  Example   : print_file($file);

=cut



sub print_file{
  my ($file) = @_;
  my $md5 = $file->md5;
  $md5 = ' ' if(!$md5);
  print join("\t", $file->name, $file->size, $md5, $file->type);
  print "\n";
}


=head2 withdraw_file

  Arg [1]   : ReseqTrack::File
  Arg [2]   : ReseqTrack::DBSQL::DBAdaptor
  Arg [3]   : string, comment for reason of withdrawal
  Arg [4]   : string, old directory
  Arg [5]   : string, new directory
  Arg [6]   : binary, if set check before running a move
  Function  : move files to new location, update those files as withdrawn in
  database
  Returntype: ReseqTrack::File 
  Exceptions: 
  Example   : 

=cut


sub withdraw_file{
  my ($file, $db, $comment, $old, $new, $check) = @_;
  $old = 'sequence_staging' unless($old);
  $new = 'sequence_staging/withdrawn' unless($new);
  my $path = $file->name;
  my $new_path = $file->name;
  $new_path =~ s/$old/$new/;
  if($path eq $new_path){
    warning("The path hasn't been changed by substitution ".$old." with ".$new.
            " have nothing to do\n");
    return undef;
  }
  my $run = 0;
  if($check){
    print "Do you want to withdraw ".$file->name." and move ".$path." to ".
        $new_path."\n";
    print "Answer Y or N\n";
    $run = get_input_arg();
  }else{
    $run = 1;
  }
  if($run){
    my $fa = $db->get_FileAdaptor;
    print "mv ".$path." ".$new_path."\n";
    move($path, $new_path);
    throw("Move of ".$path." to ".$new_path." failed") unless(-e $new_path);
    $file->withdrawn(1);
    $file->name($new_path);
    $comment = "Witdrawning ".$path." to ".$new_path unless($comment);
    my $history = ReseqTrack::History->new(
      -other_id => $old->dbID,
      -table_name => 'file',
      -comment => $comment,
        );
    $file->history($file);
    $fa->update_withdrawn($file);
  }else{
    print "You should withdrawn ".$path." to ".$new_path." and mark ".$file->dbID.
        " as withdrawn\n";
  }
  return $file;
}

sub copy_file_object{
  my ($file, $keep_md5) = @_;
  $keep_md5 = 1 unless(defined($keep_md5));
  my $new_file = ReseqTrack::File->new(
    -name => $file->name,
    -type => $file->type,
    -size => $file->size,
    -host => $file->host,
    -created => $file->created,
      );
  if($file->history && @{$file->history} >= 1){
    $new_file->history($file->history);
  }
  if($file->statistics && @{$file->statistics} >= 1){
    $new_file->statistics($file->statistics);
  }
  $new_file->md5($file->md5) if($keep_md5);
  return $new_file;
}

sub assign_type{
  my ($files) =@_;
  foreach my $file(@$files){
    if($file->name =~ /reference/i){
      $file->type("REFERENCE");
    }elsif($file->name =~ /\/technical\/working/i){
      $file->type("INTERNAL");
    }elsif($file->name =~ /changelog/i){
      $file->type("CHANGELOG");
    }elsif($file->name =~ /\/technical\//i && $file->name !~ /\/ncbi_varpipe_data\//i ){
      $file->type("TECHNICAL");
    }elsif($file->name =~ /\/release/i){
      $file->type("RELEASE");
    }elsif($file->name =~ /\.filt\.fastq\.gz$/){
      $file->type("FILTERED_FASTQ");
    }elsif($file->name =~ /sequence_staging/){
      $file->type("ARCHIVE_FASTQ");
    }elsif($file->name =~ /\.fastq\.gz/i){
      $file->type("FASTQ");
    }elsif($file->name =~ /\.bam$/){
      $file->type("BAM");
    }elsif($file->name =~ /\.bai$/){
      $file->type("BAI");
    }elsif($file->name =~ /\.bas$/ && !($file->name =~ /README/)){
      $file->type("BAS");
    }elsif($file->name =~ /\.index$/){
      $file->type("INDEX");
    }else{
      $file->type("MISC");
    }
    if($file->name =~ /\/pilot3_/){
      my $type = $file->type;
      $type = "PILOT3_".$type;
      $file->type($type);
    }elsif($file->name =~ /\/pilot_data\//){
      my $type = $file->type;
      $type = "PILOT_".$type;
      $file->type($type);
    }
    if($file->name =~ /\/withdrawn\//){
      my $type = $file->type;
      $type = "WITHDRAWN_".$type;
      $file->type($type);
    }
    if ($file->name =~ /mosaik/i ){
      my $type = $file->type;
      $type = "NCBI_".$type;
      $file->type($type);
    }
  }
  return $files;
}

sub check_type{
  my ($file) = @_;
  my $hash = standard_type_hash();
  my $type;
  return $hash->{$file->type};
}

sub standard_type_hash{
  my %hash;
  $hash{'README'} = 1;
  $hash{'CHANGELOG'} = 1;
  $hash{'INDEX'} = 1;
  $hash{'MISC'} = 1;
  $hash{'PILOT_README'} = 1;
  $hash{'RELEASE'} = 1;
  $hash{'TECHNICAL'} = 1;
  $hash{'INTERNAL'} = 1;
  $hash{'PILOT_REFERENCE'} = 1;
  $hash{'PILOT_INTERNAL'} = 1;
  $hash{'PILOT_TECHNICAL'} = 1;
  $hash{'PILOT_BAS'} = 1;
  $hash{'PILOT_BAM'} = 1;
  $hash{'PILOT_BAI'} = 1;
  $hash{'PILOT_RELEASE'} = 1;
  $hash{'PILOT_INDEX'} = 1;
  $hash{'PILOT3_BAS'} = 1;
  $hash{'PILOT3_BAM'} = 1;
  $hash{'PILOT3_BAI'} = 1;
  $hash{'REFERENCE'} = 1;
  $hash{'PILOT_INDEX'} = 1;
  $hash{'PILOT_MISC'} = 1;
  $hash{'FILTERED_FASTQ'} = 1;
  $hash{'PILOT3_FASTQ'} = 1;
  $hash{'PILOT_FASTQ'} = 1;
  $hash{'ARCHIVE_FASTQ'} = 1;
  $hash{'WITHDRAWN_FILTERED_FASTQ'} = 1;
  $hash{'WITHDRAWN_PILOT_INTERNAL'} = 1;
  $hash{'BAM'} = 1;
  $hash{'BAI'} = 1;
  $hash{'BAS'} = 1;
  $hash{'NCBI_BAM'} = 1;
  $hash{'NCBI_BAI'} = 1;
  $hash{'NCBI_BAS'} = 1;
  return \%hash;
}

=head2 move_file_in_db_and_dir

  Arg [1]   : a reference to an array of file objects to be moved
  Arg [2]   : directory name where the files should be moved to
  Arg [3]   : new file type
  Arg [4]   : db object
 
  Function  : move files to a specified directory and make the change in file table of the tracking db
  Returntype: a reference to an array of new file objects (new path and new file type); with this, new collection can be created if necessary
  Exceptions: 
  			  
  Example   : move_file_in_db_and_dir($others, $dir, $new_file_type, $db);
 
=cut

sub move_file_in_db_and_dir { 
  my ($file_objects, $dir, $new_f_type, $db) = @_;
  
  my @new_file_objects;
  
  unless ( -e $dir ) {
    mkpath($dir);
  }
  
  my $fa = $db->get_FileAdaptor;
  
  foreach my $f_obj ( @$file_objects ) {
    
    my $file_path = $f_obj->name;
    my $base_name = basename($file_path);
    
    #print "moving bam file $file_path to archive staging........\n";	
    #`mv $file_path $dir`;
    
    my $new_path = $dir . "/" . $base_name;
    move($file_path, $new_path);
    throw("Failed to move ".$file_path." to ".$new_path) unless(-e $new_path);
    $f_obj->name($new_path);
    $f_obj->type($new_f_type);
    
    my $history;
    my $hist_a = $db->get_HistoryAdaptor();
    my $possible_stored_file = $fa->fetch_by_filename($f_obj->filename); #filename function gets the basename of the file; name function gets the whole path 	
    
    if($possible_stored_file && @$possible_stored_file >0 ){
      if (@$possible_stored_file == 1) {
	#### get the possible stored object
	my $stored_file = $possible_stored_file->[0];	
	#print "matching filename stored in db: " . $stored_file->name . "\n";	
	
	#### compare new file and the exisitng file to see what the difference between the files
	my $comment = calculate_comment($f_obj, $stored_file); 
	
	$history = ReseqTrack::History->new(
					    -other_id => $stored_file->dbID,
					    -table_name => 'file',
					    -comment => $comment, 
					   );
	
	#### update the history object and assign it to the new file
	$f_obj->history($history);
	
	#### modify the file object and store it in the file table with the update function
	$f_obj->dbID($stored_file->dbID);
	$f_obj->created($stored_file->created);
	$fa->update($f_obj,1); #1 is withdraw indicator for overwrite the old file in file table 
      }
      elsif ( @$possible_stored_file > 1) {
	warning("There are multiple files in the database with the filename ".
		$f_obj->filename." Don't know what to do now");
      }
    }
    #### If the file object does not exist in db, store it and create a history object for it
    else{
      throw("Can't store ".$f_obj->dbID." ".$f_obj->name) unless(-e $f_obj->name);
      $fa->store($f_obj);	
      $history = ReseqTrack::History->new(
					  -other_id => $f_obj->dbID,
					  -table_name => 'file',
					  -comment => "Created new BAM file object", 
					 );
      $hist_a->store($history);
    }	
    push (@new_file_objects, $f_obj);					
  }#end of foreach
  
  return (\@new_file_objects);	
}	


sub get_count_stats{
  my ($file) = @_;
  my $stats = $file->statistics;
  my ($read_count, $base_count);
  foreach my $stat(@$stats){
    $read_count = $stat->attribute_value if($stat->attribute_name eq 'read_count');
    $base_count = $stat->attribute_value if($stat->attribute_name eq 'base_count');
  }
  return ($read_count, $base_count);
}

sub get_run_id_to_file_hash{
  my ($files) = @_;
  my %hash;
  foreach my $file(@$files){
    my ($run_id) = $file->filename =~ /^([E|S]RR\d+)/;
    throw("Failed to parse a run id from ".$file->name) unless($run_id);
    push(@{$hash{$run_id}}, $file);
  }
  return \%hash;
}

sub write_log {
	my ($f_obj, $log_adaptor, $message) = @_;	
	my $new_log_obj;
	
	if ($message) {
		$new_log_obj = ReseqTrack::RejectLog->new(
			-file_id 		=> $f_obj->dbID,
			-is_reject 		=> 'y',
			-reject_reason 	=> $message, 
			);
	}
	else { #a good file
		$new_log_obj = ReseqTrack::RejectLog->new(
			-file_id => $f_obj->dbID,
			-is_reject => 'n', 
			);
	}		
	
	#$log_adaptor->store($new_log_obj, 1) if ($run);
	$log_adaptor->store($new_log_obj, 1);
	
	return 1;
}	


1;
