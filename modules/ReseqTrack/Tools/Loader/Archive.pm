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
  $priority,  $from_db,

   ) = rearrange(
  [
   qw(
     LIST_FILE PATH_LIKE ACTION ACTION_LOCATION_NAME  ARCHIVE_SLEEP
     PRIORITY FROM_DB
     )
  ],
  @args
   );

 $self->list_file($list_file);
 $self->path_like($path_like);
 $self->action_string($action);
 $self->action_location_name($action_location_name);
 $self->archive_sleep($archive_sleep);
 $self->priority($priority);
 $self->from_db($from_db);
 
 $self->archive_action_adaptor();
 $self->archive_location_adaptor();
 
 $self->create_action_hash();
 $self->set_archive_actions();
 
 $self->set_location_names($action);
 $self->set_archive_location();

 $self->set_location_root();

 return $self;
}


####################################

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
####################################
sub set_archive_actions {
 my ( $self, $arg ) = @_;
 if ( !defined $self->archive_actions ) {
  my $aa              = $self->archive_action_adaptor;
  my $archive_actions = $aa->fetch_all;
  $self->archive_actions($archive_actions);
 }
}

###################################
sub archive_objects {
 my ( $self, $run ) = @_;
 my @tmp              = $self->file_paths;
 my @files_to_archive = $tmp[0];
 my %file_objects;

 my $fa = $self->db->get_FileAdaptor;
 my $aa = $self->db->get_ArchiveAdaptor;

 foreach my $file_path (@files_to_archive) {
  next unless ( -e $file_path );
  my $action = $self->{which_action_hash}{$file_path};
  $file_path =~ s/\/$//;
  my $file;

  unless ( $file_objects{$file_path} ) {
   $file = $fa->fetch_by_name($file_path);
   if ( !$file ) {
    throw( "Failed to fetch file from " . $file_path );
   }
  }
  else {
   $file = $file_objects{$file_path};
  }

  #push( @{ $changelog_hash{ $action->action } }, $file_path );
  #$changelog_name = lc( $file->type ) unless ($changelog_name);
  my $archive =
    create_archive_from_objects( $file, $action, $self->archive_location );
  $archive->priority( $self->priority );
  $aa->store($archive) if ($run);
 }

}

sub create_action_hash {
 my ( $self, $arg ) = @_;
 my %action_hash;

 my $archive_actions = $self->archive_actions;

 foreach my $action (@$archive_actions) {
  $action_hash{ $action->action } = $action;
 }
 $self->action_hash( \%action_hash );
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

#======================================

sub sanity_check_objects {
 my ( $self, $arg ) = @_;
 print "Doing Sanity check\n";


 #sanity check files
 my @files_to_archive;
 my %which_action_hash;
 my %changelog_hash;
 my @files = $self->file_paths;

 my $mini_array = $files[0];

 my $action_hash = $self->action_hash;

 foreach my $file (@$mini_array) {
  print $file, "\n";
  warning( "Can't archive a file " . $file . " which doesn't exist" )
    unless ( -e $file );

  unless ( $file =~ /$self->location_root/ ) {
   warning(   "Can't "
          . $self->action_string . " "
          . $file
          . " it isn't in the "
          . $self->location_root );
  }

  my $new_file = $file;

  my $new_root = $self->other_location->location;
  if ( $new_root =~ /\/nfs\/1000g-archive/ ) {
   $new_root .= "/vol1";
  }

  $new_file =~ s/$self->location_root/$self->new_root/;

  if ( $self->archive_location->location_name eq 'staging' ) {
   if ( -e $new_file ) {

    #print STDERR "Running replace on ".$file."\n"
    #    unless($action_string eq 'replace');
    $which_action_hash{$file} = $$action_hash{'replace'};
   }
   else {
    $which_action_hash{$file} = $$action_hash{'archive'};
   }
  }
  else {
   if ( -e $new_file ) {
    warning( "Can't dearchive " . $file . " as " . $new_file . " exists" );
    next;
   }
   else {
    $which_action_hash{$file} = $$action_hash{'dearchive'};
   }
  }
  push( @files_to_archive, $file );
 }

 $self->which_action_hash( \%which_action_hash );

 print "Done sanity check\n";

}
#######
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
#######

sub process_input {
 my ( $self, $arg ) = @_;
 my %file_objects;
 my $files;

#mystically:::  arrayref off command line
if ( $self->file ){
   $files = $self->file ;
 }


 my $fa = $self->db->get_FileAdaptor;

 if ( $self->list_file ) {
  my $list = get_lines_from_file( $self->list_file );
  push( @$files, @$list );
 }
 elsif ( $self->dir ) {
  my ( $list, $hash ) = list_files_in_dir( $self->dir, 1 );
  if ( $self->descend ) {
   push( @$files, @$list );
  }
  else {
   my $dir_list = $hash->{ $self->dir };
   push( @$files, @$dir_list );
  }
 }
 elsif ( $self->from_db ) {
  my $list;
  if ( $self->path_like && !$self->type ) {
   $list = $fa->fetch_all_like_path( $self->path_like );
  }
  elsif ( !$self->path_like && $self->type ) {
   $list = $fa->fetch_by_type( $self->type );
  }
  elsif ( $self->path_like && $self->type ) {
   my $temps = $fa->fetch_all_like_path( $self->path_like );
   foreach my $temp (@$temps) {
    next unless ( $temp->type eq $self->type );
    push( @$list, $temp );
   }
  }
  else {
   print
     "When archiving files from the database you must specify a type and/or "
     . "a root path for the files to be from\n";
  }
  
   print "Have " . @$list . " files to check\n" ;#if ( $self->verbose );

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

 $self->file_paths(@$files);
 print Dumper(@$files) if ( $self->debug );

}

############################################################

sub initial_table_clean {
 my ( $self, $arg ) = @_;
 my $db = $self->db;
 print "Doing initial table clean\n";
 my $aa       = $db->get_ArchiveAdaptor;
 my $archives = $aa->fetch_all;

 print ref($archives), "\n";
 print "Files in archive_table= ", scalar(@$archives), "\n";

 #cleanup_archive( $archives, $db, 0 );
}
##########
sub set_location_root {
 my ( $self, $arg ) = @_;
 my $archive_location = $self->archive_location;

 my $location_root = $archive_location->location;
 throw("Could not set 'location_root") if ( !$location_root );

 $self->location_root($location_root);
}

####
sub set_other_location {
 my ( $self, $arg ) = @_;

 my $other_location_adaptor = $self->archive_location_adaptor();
 my $other_location_name    = $self->archive_location_name();

 my $other_location =
   $other_location_adaptor->fetch_by_archive_location_name(
                                                          $other_location_name);

 $self->other_location($other_location);
 throw "Could not set 'other_location " if ( !$self->other_location );
}

#######

sub set_archive_location {
 my ( $self, $arg ) = @_;
 my $archive_location_adaptor = $self->archive_location_adaptor();
 my $action_location_name     = $self->action_location_name();

 my $archive_location =
   $archive_location_adaptor->fetch_by_archive_location_name(
                                                         $action_location_name);

 $self->archive_location($archive_location);

 throw "Could not set archive_location " if ( !$self->archive_location );

}
########
sub set_location_names {
 my ( $self, $arg ) = @_;

 throw("You have not set an archive action (archive replace dearchive)")
   if ( !$arg );

 if ($arg) {
  if ( $arg eq "archive" || $arg eq 'replace' ) {
   $self->action_location_name("staging");
   $self->other_location_name("archive");
  }
  elsif ( $arg eq "dearchive" ) {
   $self->action_location_name("archive");
   $self->other_location_name("staging");
  }
  else {
   throw(   "Don't know what to do with action " 
          . $arg
          . " should be archive or dearchive" );
  }
 }

 return $self->{action_string};
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
 $self->{archive_location} = $arg if ($arg);
 return $self->{archive_location};
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
}

sub archive_location {
 my ( $self, $arg ) = @_;
 $self->{archive_location} = $arg if ($arg);
 return $self->{archive_location};
}

#sub archive_action_adaptor {
# my ( $self, $arg ) = @_;
#$self->{archive_action_adaptor} = $arg if ( defined $arg );
#return $self->{archive_action_adaptor};
#}

#sub set_new_root {
# my ( $self, $arg ) = @_;

# my $other_location = $self->other_location;

# my $new_root = $other_location->location;
#if ( $new_root =~ /\/nfs\/1000g-archive/ ) {
#$new_root .= "/vol1";
#}

#throw("Could not set 'new_root") if ( !$new_root );
#$self->new_root($new_root);

#}

#sub new_root {
# my ( $self, $arg ) = @_;
#$self->{new_root} = $arg if ( defined $arg );
#return $self->{new_root};
#}

#sub create_archive_action_adaptor {
# my ( $self, $arg ) = @_;
# my $aa = "fred";

#if ( !defined $self->archive_action_adaptor || ( !defined $arg ) ) {
#print "Creating archive_action_adaptor\n";
#my $db = $self->db;
#my $aa = $db->get_ArchiveActionAdaptor;
#$self->archive_action_adaptor($aa);
#}
#elsif ( $arg->isa("ReseqTrack::ArchiveActionAdaptor") ) {
#$self->archive_action_adaptor($arg);
#}

#}
#sub create_archive_actions {
# my ( $self, $arg ) = @_;
#if ( !defined $self->archive_actions ) {
# my $aa              = $self->archive_action_adaptor;
# my $archive_actions = $aa->fetch_all;
#$self->archive_actions($archive_actions);
#}
#print Dumper ($self) if ( $self->verbose );
#}

####

#sub create_archive_location_adaptor {
# my ( $self, $arg ) = @_;
# my $la = "fred";

#if ( !defined $self->archive_location_adaptor || ( !defined $arg ) ) {
#print "Creating archive_location_adaptor\n";
#my $db = $self->db;
#$la = $db->get_ArchiveLocationAdaptor;
# $self->archive_location_adaptor($la);
#}
#elsif ( $arg->isa("ReseqTrack::ArchiveLocationAdaptor") ) {
#$self->archive_location_adaptor($arg);
#}

#print Dumper ($self) if ( $self->verbose );
#}

#sub archive_location_adaptor {
# my ( $self, $arg ) = @_;
#$self->{archive_location_adaptor} = $arg if ( defined $arg );
#return $self->{archive_location_adaptor};
#}

=pod
sub check_input {
 my $self = shift;
 
  my $action =  $self->action_string;
  
 if (! $action){
  throw "No archive action specified. Use: -action archive\/dearchive";
 }
  
   unless (     $action eq "archive" ||  $action eq "dearchive"){
     throw "Strange archive ($action) action specied. Use:archive\/dearchive";
    }
 
}
=cut

1;
