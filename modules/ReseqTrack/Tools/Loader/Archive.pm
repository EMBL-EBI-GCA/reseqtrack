package ReseqTrack::Tools::Loader::Archive;

use strict;
use warnings;

use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::ArchiveUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::GeneralUtils;
use File::Basename;
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list );
use ReseqTrack::Tools::FileUtils qw(create_history assign_type check_type);

use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(ReseqTrack::Tools::Loader);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  my (
      $list_file, $path_like, $action, $action_location_name, $archive_sleep,
      $priority,  $from_db,  $lines_check , $max_number,
     ) = rearrange(
		   [
		    qw(
     LIST_FILE PATH_LIKE ACTION ACTION_LOCATION_NAME  ARCHIVE_SLEEP
     PRIORITY FROM_DB LINES_CHECK MAX_NUMBER
     )
		   ],
		   @args
		  );
   
  $self->lines_check (1);	# default 'on';
  $self->archive_sleep ("240");

  $self->lines_check($lines_check);
  $self->list_file($list_file);
  $self->path_like($path_like);
  $self->action_string($action);
  $self->action_location_name($action_location_name);
  $self->archive_sleep($archive_sleep);
  $self->priority($priority);
  $self->from_db($from_db);
 
  $self->archive_action_adaptor();
  $self->archive_location_adaptor();


  $self->get_archive_actions(); 
  $self->create_action_hash();  
 
 
  $self->set_location_names($action); 
  $self->set_archive_location();

  $self->set_other_location();


  $self->set_location_root();
  $self->max_number ( $max_number);
 
  return $self;
}



sub cleanup_archive_table{
  my ($self, $verbose) = @_;

  my $db = $self->db;

  my $aa = $db->get_ArchiveAdaptor;

  my $archives = $aa->fetch_all;
  if (!$archives){
    print "Archive table is clean\n";
    return;
  }

  cleanup_archive( $archives, $db, $verbose) ;
  
  my $new_archives = $aa->fetch_all;
 
  return(@$new_archives);
}


sub DESTROY{
 my ( $self, $arg ) = @_;


 my $aa = $self->archive_adaptor();
 if ($aa){
   print "Deleting archive table lock\n";
   $aa->delete_archive_lock();
 }else{
   print "Failed to delete archive_lock\n";
 }

 return;

}


sub archive_adaptor{
  my ( $self, $arg ) = @_;
  if (defined $arg){
   if ( ! ($arg->isa("ReseqTrack::DBSQL::ArchiveAdaptor")) ){
     print ref ($arg), "\n";
     throw "Passed object other than ArchiveAdaptor\n"       	
   }
   $self->{archive_adaptor} = $arg
 }
 
  return $self->{archive_adaptor};
}


###################################
sub archive_objects {
  my ( $self ) = @_;
  my $tmp  = $self->file_paths;
  my $files_to_archive;
  my %file_objects;

 
  if (ref($tmp) eq 'ARRAY') {
    push(@$files_to_archive, @$tmp);
  } else {
    push(@{$files_to_archive}, $tmp);
  }

  my $total = scalar ( @{$files_to_archive});
   
  
  my $fa = $self->db->get_FileAdaptor;
  my $aa = $self->db->get_ArchiveAdaptor;
  
  print "Starting archive of $total objects\n";

  foreach my $file_path (@{$files_to_archive}) {
    print  "$file_path\n" if ($total < 10);
    
    next unless ( -e $file_path );
    my $action = $self->{which_action_hash}{$file_path};
  
    $file_path =~ s/\/$//;
    my $file;

    unless ( $file_objects{$file_path} ) {
      $file = $fa->fetch_by_name($file_path);
      if ( !$file ) {
	throw( "Failed to fetch file from " . $file_path );
      }
    } else {
      $file = $file_objects{$file_path};
    }

    my $archive =
      create_archive_from_objects( $file, $action, $self->archive_location );
    $archive->priority( $self->priority );
 
    
      has_to_many_archive_lines($self->max_number, $self->archive_sleep, $self->db) if($self->lines_check);
      $aa->store($archive);
   
   
  }

  $self->archive_adaptor($aa);

  return;
}




#======================================

sub sanity_check_objects {
  my ( $self, $arg ) = @_;
  print "Doing Sanity check\n";


  #sanity check files
  my @files_to_archive;
  my %which_action_hash;
  my %changelog_hash;
  my $files = $self->file_paths;

  my $file_array;

 
 

  if (ref($files) eq 'ARRAY') {
    push(@$file_array, @$files);
  } else {
    push(@{$file_array}, $files);
  }


  my $total_files = scalar @{$file_array};

  my $action_hash = $self->action_hash;

  foreach my $file (@$file_array) {
    warning( "Can't archive a file " . $file . " which doesn't exist" )
      unless ( -e $file );

    my $loc = $self->location_root;
    unless ( $file =~ /$loc/ ) {
      warning(   "Can't "
		 . $self->action_string . " "
		 . $file
		 . " it isn't in the "
		 . $self->location_root );
      next;
    }

    my $new_file = $file;

  
    my $new_root = $self->other_location->location;
    if ( $new_root =~ /\/nfs\/1000g-archive/ ) {
      $new_root .= "/vol1";
    }

    my $loc_root = $self->location_root;
    $new_file =~ s/$loc_root/$new_root/;

    

    if ( $self->archive_location->location_name eq 'staging' ) {

      if ( -e $new_file ) {
	print STDERR "Running replace on: ".$file."\n" if ($self->verbose);
	#    unless($action_string eq 'replace');
	$which_action_hash{$file} = $$action_hash{'replace'};
	
      } else {
	print STDERR "Running archive on: ".$file."\n" if ($self->verbose);
	$which_action_hash{$file} = $$action_hash{'archive'};
      }
    } else {
      if ( -e $new_file ) {
	warning( "Can't dearchive " . $file . " as " . $new_file . " exists" );
	next;
      } else {
	$which_action_hash{$file} = $$action_hash{'dearchive'};
      }
    }
    push( @files_to_archive, $file );
  }

  $self->which_action_hash( \%which_action_hash );

  print "Done sanity check\n" if $self->verbose;
  return;
}
#######


sub process_input {
  my ( $self, $arg ) = @_;
  my %file_objects;
  my $files;

  #mystically:::  arrayref off command line
  if ( $self->file ) {
    $files = $self->file ;
  }


  my $fa = $self->db->get_FileAdaptor;

  if ( $self->list_file ) {
    my $list = get_lines_from_file( $self->list_file );
    push( @$files, @$list );
  } elsif ( $self->dir ) {
    my ( $list, $hash ) = list_files_in_dir( $self->dir, 1 );
    if ( $self->descend ) {
      push( @$files, @$list );
    } else {
      my $dir_list = $hash->{ $self->dir };
      push( @$files, @$dir_list );
    }
  } elsif ( $self->from_db ) {
    my $list;
    if ( $self->path_like && !$self->type ) {
      $list = $fa->fetch_all_like_path( $self->path_like );
    } elsif ( !$self->path_like && $self->type ) {
      $list = $fa->fetch_by_type( $self->type );
    } elsif ( $self->path_like && $self->type ) {
      my $temps = $fa->fetch_all_like_path( $self->path_like );
      foreach my $temp (@$temps) {
	next unless ( $temp->type eq $self->type );
	push( @$list, $temp );
      }
    } else {
      print
	"When archiving files from the database you must specify a type and/or "
	  . "a root path for the files to be from\n";
    }
  
    print "Have " . @$list . " files to check\n" ; #if ( $self->verbose );

    foreach my $file (@$list) {
      my $other_root = $self->other_location->location;
      next if ( $file->path =~ /$other_root/ );
      push( @$files, $file->full_path );
      $file_objects{ $file->full_path } = $file;
    }
    print "Have " . @$list . " file objects  to process\n";
  }

  print "Have " . @$files . " file paths  to process\n";

  throw(   "Need more than zero froms in file array from either the -file, "
	   . "-file_list or -dir" )
    unless ( @$files >= 1 );

  $self->file_paths($files);
  print Dumper($files) if ( $self->debug );
  return;
}

############################################################
sub which_action_hash {
  my ( $self, $arg ) = @_;
  if ( defined $arg && ref($arg) ne "HASH" ) {
    throw "Bad pass. to action_hash. Not passing HASH ref";
  }
  if ( ref($arg) eq "HASH" && ( defined $arg ) ) {
    $self->{which_action_hash} = $arg;
  }
  return $self->{which_action_hash};
}
sub set_other_location {
  my ( $self, $arg ) = @_;

  my $ala = $self->archive_location_adaptor();
  
  my $archive_locations = $ala->fetch_all;

  foreach my $location (@$archive_locations) {
    if ( $location->location_name eq $self->other_location_name) {
      $self->other_location ($location);
      return;
    }
  }

  throw "Could not set other_location " if ( !$self->other_location );
  return;
}
####################################
sub get_archive_actions {
  my ( $self, $arg ) = @_;
  if ( !defined $self->archive_actions ) {
    my $aa              = $self->archive_action_adaptor;
    my $archive_actions = $aa->fetch_all;
    $self->archive_actions($archive_actions);
  }
  return;
}
####################################
sub create_action_hash {
  my ( $self, $arg ) = @_;
  my %action_hash;

  my $archive_actions = $self->archive_actions;

  foreach my $action (@$archive_actions) {
    $action_hash{ $action->action } = $action;
  }
  $self->action_hash( \%action_hash );
  return;
}
sub archive_location_adaptor {
  my ( $self, $arg ) = @_;

  if ($arg) {
    throw(
	  "Must pass archive_action_adaptor a 
          ReseqTrack::DBSQL::ArchiveLocationAdaptor not" . $arg
	 ) unless ( $arg->isa("ReseqTrack::DBSQL::ArchiveLocationAdaptor") );

    $self->{archive_location_adaptor} = $self->db->get_ArchiveLocationAdaptor;
  }

  unless ( $self->{archive_location_adaptor} ) { 
    $self->{archive_location_adaptor} = $self->db->get_ArchiveLocationAdaptor;
  }

  return $self->{archive_location_adaptor};
}
####################################
sub archive_action_adaptor {
  my ( $self, $arg ) = @_;

  if ($arg) {
    throw(
	  "Must pass archive_action_adaptor a 
          ReseqTrack::DBSQL::ArchiveActionAdaptor not" . $arg
	 ) unless ( $arg->isa("ReseqTrack::DBSQL::ArchiveActionAdaptor") );

    $self->{archive_action_adaptor} = $arg;
  }
  unless ( $self->{archive_action_adaptor} ) {
    $self->{archive_action_adaptor} = $self->db->get_ArchiveActionAdaptor;
  }

  return $self->{archive_action_adaptor};
}

###############################
sub set_location_names {
  my ( $self, $arg ) = @_;

  throw("You have not set an archive action (archive replace dearchive)")
    if ( !$arg );

  if ($arg) {
    if ( $arg eq "archive" || $arg eq 'replace' ) {
      $self->action_location_name("staging");
      $self->other_location_name("archive");
    } elsif ( $arg eq "dearchive" ) {
      $self->action_location_name("archive");
      $self->other_location_name("staging");
    } else {
      throw(   "Don't know what to do with action " 
	       . $arg
	       . " should be archive or dearchive" );
    }
 
  }
 
  return $self->{action_string};
}
#######

sub set_archive_location {
  my ( $self, $arg ) = @_;
  my $ala = $self->archive_location_adaptor();
  
  my $archive_locations = $ala->fetch_all;

  foreach my $location (@$archive_locations) {
    if ( $location->location_name eq $self->action_location_name) {
      $self->archive_location ($location);
      return;
    }
  }

  throw "Could not set archive_location " if ( !$self->archive_location );
  return;
}

sub action_hash {
  my ( $self, $arg ) = @_;
  if ( defined $arg && ref($arg) ne "HASH" ) {
    throw "Bad pass. to action_hash. Not passing HASH ref";
  }

  if ( ref($arg) eq "HASH" && ( defined $arg ) ) {
    $self->{action_hash} = $arg;
  }
  return $self->{action_hash};
}


##########
sub set_location_root {
  my ( $self, $arg ) = @_;
  my $archive_location = $self->archive_location;

  my $location_root = $archive_location->location;
  throw("Could not set 'location_root") if ( !$location_root );

  $self->location_root($location_root);
  return;
}



sub archive_actions {
  my ( $self, $arg ) = @_;
  $self->{archive_actions} = $arg if ($arg);
  return $self->{archive_actions};
}

sub location_root {
  my ( $self, $arg ) = @_;
  $self->{location_root} = $arg if ($arg);
  return $self->{location_root};
}

sub path_like {
  my ( $self, $arg ) = @_;
  $self->{path_like} = $arg if ($arg);
  return $self->{path_like};
}

sub action_string {
  my ( $self, $arg ) = @_;
  $self->{action_string} = $arg if ($arg);
  return $self->{action_string};
}

sub action_location_name {
  my ( $self, $arg ) = @_;
  $self->{action_location_name} = $arg if ($arg);
  return $self->{action_location_name};
}

sub other_location_name {
  my ( $self, $arg ) = @_;
  $self->{other_location_name} = $arg if ($arg);
  return $self->{other_location_name};
}

sub other_location {
  my ( $self, $arg ) = @_;
  $self->{other_location} = $arg if ($arg);
  return $self->{other_location};
}

sub archive_sleep {
  my ( $self, $arg ) = @_;
  $self->{sleep} = $arg if ($arg);
  return $self->{sleep};
}

sub from_db {
  my ( $self, $arg ) = @_;
  $self->{from_db} = $arg if ($arg);
  return $self->{from_db};
}

sub priority {
  my ( $self, $arg ) = @_;
  $self->{priority} = $arg if ($arg);
  return $self->{priority};
}

sub archive_location {
  my ( $self, $arg ) = @_;
  $self->{archive_location} = $arg if ($arg);
  return $self->{archive_location};
}


sub lines_check{
  my ( $self, $arg ) = @_;
  $self->{lines_check} = $arg if (defined $arg);
  return $self->{lines_check};
}

sub max_number{
  my ( $self, $arg ) = @_;
  $self->{max_number} = $arg if (defined $arg);
  return $self->{max_number};
}



1;
