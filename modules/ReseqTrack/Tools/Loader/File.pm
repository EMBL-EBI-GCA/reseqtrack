package ReseqTrack::Tools::Loader::File;
use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list );
use ReseqTrack::Tools::FileUtils qw(create_history assign_type check_type);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::DBSQL::DBAdaptor;
use Data::Dumper;



sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
 
 
 my (
      $db,          $dbhost,           $dbname,
      $dbuser,      $dbpass,           $dbport,
      $dir,         $use_dir,          $file,
      $list_file,   $md5_file,         $type,
      $host,        $run,              $verbose,
      $descend,     $die_for_problems, $update_existing,
      $store_new,   $help,             $assign_types,
      $check_types, $input,
   )
   = rearrange(
  [
   qw(
     DB     DBHOST DBNAME DBUSER
     DBPASS DBPORT DIR
     USE_DIR FILE LIST_FILE
     MD5_FILE
     TYPE HOST RUN
     VERBOSE DESCEND DIE_FOR_PROBLEMS
     UPDATE_EXISTING STORE_NEW
     HELP  ASSIGN_TYPES
     CHECK_TYPES INPUT)
  ],
  @args
   );
 $self->dbhost($dbhost);
 $self->dbname($dbname);
 $self->dbuser($dbuser);
 $self->dbpass($dbpass);
 $self->dbport($dbport);
 $self->dir($dir);
 $self->use_dir($use_dir);
 $self->file($file);
 $self->list_file($list_file);
 $self->md5_file($md5_file);
 $self->type($type);
 $self->host($host);
 $self->run($run);
 $self->verbose($verbose);
 $self->descend($descend);
 $self->die_for_problems($die_for_problems);
 $self->update_existing($update_existing);
 $self->store_new($store_new);

 $self->help($help);
 $self->assign_types($assign_types);
 $self->check_types($check_types);

 $self->create_db_adaptor;    #ReseqTrack::DBSQL::DBAdaptor
                              #$self->assign_host_object;

 return $self;
}

###############################################################
sub process_input {

 my $self    = shift;
 my $use_dir = $self->use_dir;
 my $dir     = $self->dir;
 my @paths
 ;
 if ( scalar( @{ $self->file } ) ) {
  print "Load files from command line\n";
  push( @paths, $self->add_files_from_var );
 }
 elsif ( $self->list_file ) {
  print "Loading from list\n";
  push( @paths, $self->add_files_from_list_file )    # -list_file
 }

 elsif ( $self->dir ) {
  print "Load from dir\n";
  push( @paths, $self->add_files_from_dir );         # -dir
 }

 elsif ( $self->md5_file ) {
  push( @paths, $self->add_full_paths_to_md5 );
 }    # -md5 file

 throw
"Found 0 files specifed using standard options: -file -dir -md5_file -list_file"
   unless (@paths);

 $self->file_paths( \@paths );

 print "Have " . @paths . " file paths to sanity check\n";
 $self->file_sanity_check;

 print "Have " . @paths . " file paths to load\n";

}
#################################################################
sub create_file_obj_array {
 my $self  = shift;
 my $paths = $self->file_paths;
 my $type  = $self->type;
 my $host  = $self->host;

 my $files = create_objects_from_path_list( $paths, $type, $host );

 print Dumper ($files);
}

#################################################################
sub assign_host_object {
 my $self      = shift;
 my $host_name = $self->host;
 my $db        = $self->db;

 my $host = get_host_object( $host_name, $db );

 if ( !$host ) {
  throw "Failed to create host object";
 }

 $self->host_obj = $host;
 print "Have host object\n";
 print Dumper ($self);
}

sub create_db_adaptor {
 my $self = shift;
 my $db   = $self->db;

 if ( !$db ) {
  $db =
    ReseqTrack::DBSQL::DBAdaptor->new(
                                       -host   => $self->dbhost,
                                       -user   => $self->dbuser,
                                       -port   => $self->dbport,
                                       -dbname => $self->dbname,
                                       -pass   => $self->dbpass,
    );
 }

 $self->db($db);
 throw("Failed to create db adaptor") unless $self->db;
print "Have db adaptor\n";
}

#retrieve entries and correct if necessary for full paths;
sub add_full_paths_to_md5 {
 my $self    = shift;
 my $use_dir = $self->use_dir;
 my $dir     = $self->dir;
 my $file    = $self->md5_file;
 my @paths;

 my $hash = get_md5hash($file);

 foreach my $key ( keys(%$hash) ) {
  #print $key, "\n";
  my $md5       = $hash->{$key};
  my $full_path = $key;
  if ( $use_dir && $dir ) {
   $full_path = $dir . "/" . $key;
   $full_path =~ s/\/\//\//g;
  }
  $hash->{$full_path} = $md5;
  push( @paths, $full_path );
 }

 $self->md5_hash($hash);

 return @paths;
}

#########################
######################################################
######################################################

#  extract to base class

# if -dir option is involked .... get files in that dir.
sub add_files_from_dir {
 my $self    = shift;
 my $dir     = $self->dir;
 my $descend = $self->descend;
 my @paths;

 return if ( !$dir );

 #get files from dir
 if ($dir) {
  my ( $paths, $hash ) = list_files_in_dir( $dir, 1 );
  unless ($descend) {
   $paths = $hash->{$dir};
  }
  warning( "No paths were found in " . $dir ) if ( !$paths || @$paths == 0 );
  push( @paths, @$paths ) if ( $paths && @$paths >= 1 );
 }

 return @paths;

}

sub add_files_from_list_file {
 my $self      = shift;
 my $list_file = $self->list_file;
 my $use_dir   = $self->use_dir;
 my $dir       = $self->dir;
 my $full_path;
 my @paths;

 return if ( !$list_file );

 my $lines = get_lines_from_file($list_file);

 foreach my $line (@$lines) {
  my $full_path = $line;
  if ( $use_dir && $dir ) {
   $full_path = $dir . "/" . $line;
   $full_path =~ s/\/\//\//g;
  }
  next if ( -d $full_path );
  push( @paths, $full_path );
 }

 return @paths;

}

sub add_files_from_var {
 my $self    = shift;
 my $files   = $self->file;
 my $use_dir = $self->use_dir;
 my $dir     = $self->dir;
 my $full_path;
 my @paths;

 return if ( !$files );

 foreach my $file (@$files) {

  #print $file,"\n";
  my $full_path = $file;
  if ( $use_dir && $dir ) {
   $full_path = $dir . "/" . $file;
   $full_path =~ s/\/\//\//g;

   #print $full_path,"\n";

  }
  push( @paths, $full_path );
 }

 return @paths;
}

#################################################################

sub file_sanity_check {
 my $self = shift;

 my $files = $self->file_paths;

 throw "No files specifed" unless $files;

 foreach my $f (@$files) {
  my @a = split /\//, $f;

  #should never happen ..  if you have any '.' directories
  foreach (@a) {
   print("Have \. type file or dir entry\:$f\n") if ( $_ =~ /^\./ );
  }

  # better be there
  if ( !-e $f ) {
   throw "$f does not exist\n";
   next;
  }
  throw "$f does not have full path\n" unless ( $f =~ /^\// );

  #better have size > 0
  my $size = -s $f;
  print "$f has 0 size \n" unless ( $size > 0 );
 }

 #better be unique.
 my %unique;
 my %dups;
 my $tot_dups =0;
 
 foreach my $f (@$files) {
 
   if (defined $unique{$f}){
      $tot_dups++;
    $dups{$f} = 1;
    next;
     }
     $unique{$f} = 1;
 }



if ($tot_dups){
 print "Duplicate file names:\n";
 foreach my $i (keys %dups){
  print $i,"\n";  
  }
throw "Have $tot_dups duplicate file names.Please correct";
}
 
 
}

sub file_paths {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{file_paths} = $arg;
 }
 return $self->{file_paths};
}

sub db {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{db} = $arg;
 }
 return $self->{db};
}

sub md5_hash {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{md5_hash} = $arg;
 }
 return $self->{md5_hash};
}

sub check_types {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{check_types} = $arg;
 }
 return $self->{check_types};
}

sub dbhost {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{dbhost} = $arg;
 }
 return $self->{dbhost};
}

sub dbname {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{dbname} = $arg;
 }
 return $self->{dbname};
}

sub dbuser {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{dbuser} = $arg;
 }
 return $self->{dbuser};
}

sub dbpass {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{dbpass} = $arg;
 }
 return $self->{dbpass};
}

sub dbport {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{dbport} = $arg;
 }
 return $self->{dbport};
}

sub dir {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{dir} = $arg;
 }
 return $self->{dir};
}

sub use_dir {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{use_dir} = $arg;
 }
 return $self->{use_dir};
}

sub file {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{file} = $arg;
 }
 return $self->{file};
}

sub list_file {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{list_file} = $arg;
 }
 return $self->{list_file};
}

sub md5_file {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{md5_file} = $arg;
 }

 return $self->{md5_file};
}

sub type {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{type} = $arg;
 }
 return $self->{type};
}

sub host {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{host} = $arg;
 }
 return $self->{host};
}

sub run {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{run} = $arg;
 }
 return $self->{run};
}

sub verbose {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{verbose} = $arg;
 }
 return $self->{verbose};
}

sub descend {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{descend} = $arg;
 }
 return $self->{descend};
}

sub die_for_problems {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{die_for_problems} = $arg;
 }
 return $self->{die_for_problems};
}

sub update_existing {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{update_existing} = $arg;
 }
 return $self->{update_existing};
}

sub store_new {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{store_new} = $arg;
 }
 return $self->{store_new};
}

sub help {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{help} = $arg;
 }
 return $self->{help};
}

sub assign_types {
 my ( $self, $arg ) = @_;
 if ($arg) {
  $self->{assign_types} = $arg;
 }
 return $self->{assign_types};
}
1;

