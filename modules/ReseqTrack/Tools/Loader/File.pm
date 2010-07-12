package ReseqTrack::Tools::Loader::File;
use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list );
use ReseqTrack::Tools::FileUtils qw(create_history assign_type check_type );
use ReseqTrack::Tools::FileSystemUtils
  qw( get_lines_from_file run_md5 get_md5hash);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::DBSQL::DBAdaptor;
use Data::Dumper;
use File::Basename;

use vars qw(@ISA);
@ISA = qw(ReseqTrack::Tools::Loader);

sub new {
 my ( $class, @args ) = @_;
 my $self = $class->SUPER::new(@args);

 my (
      $md5_file,         $type,            $hostname,
      $die_for_problems, $update_existing, $store_new,
      $assign_types,     $check_types,     $check_md5,
      $md5_program,      $remote,
   )
   = rearrange(
  [
   qw(
     MD5_FILE     TYPE
     HOSTNAME     DIE_FOR_PROBLEMS    UPDATE_EXISTING
     STORE_NEW      ASSIGN_TYPES
     CHECK_TYPES     CHECK_MD5
     MD5_PROGRAM REMOTE
     )
  ],
  @args
   );

 #Defaults
 $self->assign_types('1')     unless ( defined $assign_types );
 $self->check_types('1')      unless ( defined $check_types );
 $self->md5_program("md5sum") unless ( defined $md5_program );

 $self->assign_types($assign_types) if ( !defined $self->assign_types );
 $self->check_types($check_types)   if ( !defined $self->check_types );
 $self->md5_program("md5sum")       if ( !defined $self->md5_program );
 $self->md5_file($md5_file);
 $self->type($type) if ( !defined $self->type );
 $self->hostname($hostname);
 $self->die_for_problems($die_for_problems);
 $self->update_existing($update_existing);
 $self->store_new($store_new);
 $self->remote($remote);
 $self->check_md5($check_md5);

 $self->assign_host_object();
 print Dumper($self) if $self->debug;
 return $self;
}

###############################################################################
sub process_input {

 my $self = shift;
 my @paths;

 if ( !$self->isa("ReseqTrack::Tools::Loader") ) {
  throw "Must pass: process_input : method a \"ReseqTrack::Tools::Loader\"";
 }

 my $inputs = 0;

 #only allowing one input mode at the moment.
 $inputs++ if ( scalar( @{ $self->file } ) );
 $inputs++ if ( $self->list_file );
 $inputs++ if ( $self->dir );
 $inputs++ if ( $self->md5_file );

=pod
 if ( $inputs > 1 ) {
  print '=' x 45;
  print "\nYou are using:\n";
  print " -file\n"      if ( scalar( @{ $self->file } > 0 ) );
  print " -list_file\n"    if ( $self->list_file );
  print " -dir\n"           if ( $self->dir );
  print " -md5_file\n"  if ( $self->md5_file );
  throw("Only allowing one type of input mode. You have $inputs");
 }
=cut

 if ( $inputs == 1 && $self->md5_file ) {
  $self->load_md5_only('1');
 }
 else {
  $self->load_md5_only('0');
 }

 if ( scalar( @{ $self->file } ) ) {
  print "Load files from command line\n";
  push( @paths, $self->add_files_from_cmd_line );
 }
 if ( $self->list_file ) {
  print "Loading from list\n";
  push( @paths, $self->add_files_from_list_file )    # -list_file
 }
 if ( $self->dir ) {
  print "Load from dir\n";
  push( @paths, $self->add_files_from_dir );         # -dir
 }
 if ( $self->md5_file ) {
  push( @paths, $self->load_md5_file );              # -md5 file
 }

 throw
"Found 0 files specifed using standard options: -file -dir -md5_file -list_file"
   unless (@paths);

 $self->file_paths( \@paths );

 print "Have " . @paths . " file paths to sanity check\n";
 $self->file_sanity_check;
 print "Have " . @paths . " files to load\n";

}
####
sub load_md5_file {
 my $self = shift;
 my @paths;
 my $path_warnings = 0;
 my $new_hash      = {};

 my $hash = get_md5hash( $self->md5_file );

 #load from md5 file only .
 if ( $self->load_md5_only ) {
  foreach my $key ( keys(%$hash) ) {
   my $md5       = $hash->{$key};
   my $full_path = $key;
   $hash->{$full_path} = $md5;
   push( @paths, $full_path );
  }

  $self->md5_hash($hash);
  return @paths;
 }

 # check md5 from file, but give file list elsewhere.
 # use only basename as reference. No duplicates in hash.
 if ( !$self->load_md5_only ) {

  foreach my $key ( keys(%$hash) ) {
   my $md5  = $hash->{$key};
   my $path = $key;
   if ( !-e $path ) {
    print "Using basename for file $key\n";
    $path = basename($key);
   }
   if ( !exists $new_hash->{$path} ) {
    $new_hash->{$path} = $md5;
   }
   else {
    throw
      "Got duplicate basenames $key in md5_file for files that do not exist ";
   }
  }
  $self->md5_hash($new_hash);
  return;
 }

}

sub create_files {
 my $self = shift;

 throw "No type specified ( -type ) and assign_types = 0"
   if ( $self->type eq "MUST_FIX" && $self->assign_types == 0 );

 my $file_objs =
   create_objects_from_path_list( $self->file_paths, $self->type, $self->host );

 if ( $self->assign_types ) {
  print "Assigning types\n";
  $file_objs = assign_type($file_objs);
 }
 else {
  print "Skipping assign_types\n";
 }

 if ( $self->check_types ) {
  print "Checking types\n";
  my @wrong;
  foreach my $file (@$file_objs) {
   push( @wrong, $file ) unless ( check_type($file) );
  }

  print STDERR "There are " . @wrong . " files with the wrong type\n"
    if ( $self->verbose );

  foreach my $file (@wrong) {
   print $file->name . " " . $file->type . " is wrong\n";
  }
  throw("There are problems with the file types") if ( @wrong >= 1 );
 }

 my $objs = scalar(@$file_objs);
 print "Created $objs file objects\n";

 $self->files($file_objs);

}

sub load_files {
 my ( $self, $run ) = @_;
 my $md5_hash = $self->md5_hash;
 my $files    = $self->files;

 print "Starting Load process\n";
 print "Not running ... really\n" unless ($run);

 #sanity checks not sure if there should be any but check for existance if
 #remote isn't specified
 my @problems;
 foreach my $file (@$files) {

  #print $file->full_path,"\n";
  my $string = $file->full_path . " doesn't exist"
    unless ( $self->remote || -e $file->full_path );
  if ( !$self->remote && -d $file->full_path ) {
   $string = $file->full_path . " is a directory ";
  }
  push( @problems, $string ) if ($string);
}

 if (@problems) {
   foreach my $problem (@problems) {
     print STDERR $problem;
   }
   if ( $self->die_for_problems ) {
     throw( @problems . " problems identified with input set dying" );
   }
 }

 if ( $self->check_md5 ) {
  warning(
   "You are going to run md5s for " . @$files . " files this may take a while" )
    if ( @$files >= 50 );
 }

 my $fa = $self->db->get_FileAdaptor;

 $md5_hash = $self->md5_hash;
 

 my @storage_problems;
FILE: foreach my $file (@$files) {
  print "Loading:", $file->full_path, "\n" if ($self->verbose);

  if ( $self->check_md5 ) {
   my $md5 = run_md5( $file->full_path, $self->md5_program );    # if ($run);
   $file->md5($md5);
  }

  if ( $md5_hash && keys %$md5_hash ) {
   my $md5 = $md5_hash->{ $file->full_path };

   if ( !$md5 ) {
    my $bname = basename( $file->full_path );    # from foreign md5 file
    $md5 = $md5_hash->{$bname};
   }
   if ( $file->md5 ) {
    unless ( $file->md5 eq $md5 ) {
     print $file->full_path
       . " has a different md5 to the once specified in "
       . "input Skipping the file\n";
     next FILE;
    }
   }
   $file->md5($md5);
  }
  eval {
   if ( $self->update_existing )
   {
    my $existing = $fa->fetch_by_name( $file->name );
    if ($existing) {
     $file->dbID( $existing->dbID );
     my $history = create_history( $file, $existing );
     $file->history($history) if ($history);
     unless ($history) {
      next FILE;
     }
    }
    else {
     my $possible_existing = $fa->fetch_by_filename( $file->filename );
     if ($possible_existing) {
      if ( @$possible_existing == 1 ) {
       my $existing = $possible_existing->[0];
       $file->dbID( $existing->dbID );
       my $history = create_history( $file, $existing );
       next FILE unless ($history);
       $file->history($history) if ($history);
      }
      elsif ( @$possible_existing >= 2 ) {
       my $for_update;
       foreach my $existing (@$possible_existing) {
        if ( $existing->type eq $file->type ) {
         if ($for_update) {
          warning(   "Can't update "
                   . $file->filename
                   . " there are multiple files "
                   . "which share its name and type" );
         }
         $for_update = $existing;
        }
       }
       $file->dbID( $for_update->dbID ) if ($for_update);
       my $history = create_history( $file, $for_update ) if ($for_update);
       next FILE unless ($history);
       $file->history($history) if ($history);
       unless ($for_update) {
        print STDERR "There are "
          . @$possible_existing
          . " possible existing files\n";
        foreach my $file (@$possible_existing) {
         print STDERR $file->dbID . " " . $file->name . "\n";
        }
        throw(   "Have multiple files linked to "
               . $file->filename
               . " not sure how "
               . "to update the file" );
       }
      }
     }
    }
   }
   $fa->store( $file, $self->update_existing, $self->store_new ) if ($run);
  };
  if ($@) {
   throw( "Problem storing " . $file . " " . $file->full_path . " $@" )
     if ( $self->die_for_problems );
   push( @storage_problems, $file->full_path . " " . $@ );
  }
 }

 foreach my $problem (@storage_problems) {
  print STDERR $problem . "\n";
 }

}
####
sub die_for_problems {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{die_for_problems} = $arg;
 }
 return $self->{die_for_problems};
}

sub update_existing {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{update_existing} = $arg;
 }
 return $self->{update_existing};
}

sub store_new {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{store_new} = $arg;
 }
 return $self->{store_new};
}

sub check_md5 {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{check_md5} = $arg;
 }
 return $self->{check_md5};
}

sub load_md5_only {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{load_md5_only} = $arg;
 }
 return $self->{load_md5_only};
}

sub md5_program {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{md5_program} = $arg;
 }
 return $self->{md5_program};
}

sub assign_types {
 my ( $self, $arg ) = @_;

 if ( defined $arg ) {
  $self->{assign_types} = $arg;
 }

 return $self->{assign_types};
}

sub check_types {
 my ( $self, $arg ) = @_;

 if ( defined $arg ) {
  $self->{check_types} = $arg;
 }
 return $self->{check_types};
}

sub remote {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{remote} = $arg;
 }
 return $self->{remote};
}

sub assign_host_object {
 my $self = shift;

 my $host = get_host_object( $self->hostname, $self->db );

 if ( !$host ) {
  throw "Failed to create host object";
 }

 $self->host($host);
 print Dumper ( $self->host ) if ( $self->debug );

}

sub file_sanity_check {
 my $self  = shift;
 my $files = $self->file_paths;
 my %bad;

 throw "No files specifed" unless $files;

 foreach my $f (@$files) {
  my @a = split /\//, $f;

  #should never happen ..  if you have any '.' directories
  foreach my $j (@a) {

   #print $j,"\n";
   throw "Have \. type file: $f    $_" if ( $j =~ /^\./ );
  }

  if ( !-e $f ) {
   $bad{$f} = "Does not exist    :";
   print "Does not exist:$f\n";
   next;
  }

  if ( -d $f ) {
   $bad{$f} = "Is directory        :";
   print "Is directory:$f\n";
   next;
  }
  
  if ( !   ($f  =~  /^\//)   ) {
   print "Not full path $f\n";
   $bad{$f} = "Full path bad      :";
   next;
  }
  #my $size = -s $f;
  #$bad{$f} = "File has 0 size:" unless ( $size > 0 );
 }

 #better be unique.
 my %unique;
 my %dups;
 my $tot_dups = 0;

 foreach my $f (@$files) {
  $dups{$f}++;
 }

 foreach my $i ( keys %dups ) {
  if ( $dups{$i} > 1 ) {
   $bad{$i} = "Duplicate file name:";
  }
 }
 if ( keys %bad ) {
  print STDERR "Found the following problems:\n";
  foreach my $i ( keys %bad ) {
   print STDERR $bad{$i}, $i, "\n";
  }
  throw "Fix file path problems";
 }
}

sub files {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{file_objs} = $arg;
 }
 return $self->{file_objs};
}

sub host {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{host} = $arg;
 }
 return $self->{host};
}
###
sub hostname {
 my ( $self, $arg ) = @_;
 if ( defined $arg ) {
  $self->{hostname} = $arg;
 }
 return $self->{hostname};
}
###
sub md5_file {
 my ( $self, $arg ) = @_;
 if (defined $arg) {
  $self->{md5_file} = $arg;
 }
 return $self->{md5_file};
}
###
sub md5_hash {
 my ( $self, $arg ) = @_;
 if (defined $arg) {
  $self->{md5_hash} = $arg;
 }
 return $self->{md5_hash};
}
1;

