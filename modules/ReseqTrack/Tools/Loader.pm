
=head1 NAME

ReseqTrack::Tools::Loader

=head1 SYNOPSIS

Misguided attempt to get Loader stuff to work
=head1 Example


=cut

package ReseqTrack::Tools::Loader;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileUtils
  qw(create_objects_from_path_list create_history assign_type check_type);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw( get_lines_from_file get_md5hash list_files_in_dir);
use Data::Dumper;

sub new {
 my ( $class, @args ) = @_;
 my $self = {};
 bless $self, $class;

 my (
      $db,     $dbhost,  $dbname,  $dbuser, $dbpass,
      $dbport, $dir,     $file,    $type, $list_file,
      $verbose, $descend,   $debug,
   )
   = rearrange(
  [
   qw(
     DB     DBHOST     DBNAME     DBUSER     DBPASS
     DBPORT     DIR     FILE     TYPE  LIST_FILE
     VERBOSE     DESCEND     DEBUG)
  ],
  @args
   );



 print "Inheriting from Loader\n" if $debug;
  #defaults
  $self->type("MUST_FIX")              unless( defined $type  ) ;
   $self->descend('1')                      unless (defined $descend);
  
   
 $self->dbhost($dbhost);
 $self->dbname($dbname);
 $self->dbuser($dbuser);
 $self->dbpass($dbpass);
 $self->dbport($dbport);
 $self->dir($dir);
 $self->file($file);
 $self->type($type);
 $self->list_file($list_file);
 $self->verbose($verbose);
 $self->descend($descend);
 $self->debug($debug);
 

 
 $self->create_DBAdaptor;    #ReseqTrack::DBSQL::DBAdaptor

 return $self;
}
###
sub process_input {

 my $self    = shift;
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

 if ( $inputs > 1 ) {
  print '=' x 45;
  print "\nYou are using:\n";
  print " -file\n"          if ( scalar( @{ $self->file } > 0 ) );
  print " -list_file\n"    if ( $self->list_file );
  print " -dir\n"           if ( $self->dir );
  print " -md5_file\n"  if ( $self->md5_file );
  throw("Only allowing one type of input mode. You have $inputs");
 }

 if ( scalar( @{ $self->file } ) ) {
  print "Load files from command line\n";
  push( @paths, $self->add_files_from_cmd_line );
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
  push( @paths, $self->add_full_paths_to_md5 ); # -md5 file
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
sub create_DBAdaptor {
 my $self = shift;
 my $db;

 if (  ! $self->db) {
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
 print Dumper ( $self->db($db) ) if $self->debug;
}
###
# if -dir option is involked .... get files in that dir.
sub add_files_from_dir {
 my $self    = shift;
 my $descend = $self->descend;
 my @paths;

 return if ( !$self->dir );
###
 #get files from dir
 if ($self->dir) {
  my ( $paths, $hash ) = list_files_in_dir( $self->dir, 1 );
  unless ($descend) {
   $paths = $hash->{$self->dir};
  }
  warning( "No paths were found in " . $self->dir ) if ( !$paths || @$paths == 0 );
  push( @paths, @$paths ) if ( $paths && @$paths >= 1 );
 }
 return @paths;
}
###
sub add_files_from_list_file {
 my $self      = shift;
  my $full_path;
  my @paths;

 return if ( !$self->list_file );

 my $lines = get_lines_from_file($self->list_file);

 foreach my $line (@$lines) {
  my $full_path = $line;
   next if ( -d $full_path );
  push( @paths, $full_path );
 }
 return @paths;
}
###
sub add_files_from_cmd_line {
 my $self    = shift;
 my $full_path;
 my @paths;

 return if ( !$self->file);

 foreach my $file (@{$self->file} ) {
  my $full_path = $file;
   push( @paths, $full_path );
 }
 return @paths;
}
####

sub file_paths {
 my ( $self, $arg ) = @_;
  $self->{file_paths} = $arg  if (defined $arg);
 return $self->{file_paths};
}
###
sub db {
 my ( $self, $arg ) = @_;
  $self->{db} = $arg if (defined $arg) ;
 return $self->{db};
}
###
sub dbhost {
 my ( $self, $arg ) = @_;
   $self->{dbhost} = $arg if (defined $arg);
 return $self->{dbhost};
}
###
sub dbname {
 my ( $self, $arg ) = @_;
   $self->{dbname} = $arg if (defined $arg);
  return $self->{dbname};
}
###
sub dbuser {
 my ( $self, $arg ) = @_;
   $self->{dbuser} = $arg if (defined $arg);
  return $self->{dbuser};
}
###
sub dbpass {
 my ( $self, $arg ) = @_;
  $self->{dbpass} = $arg if (defined $arg);
 return $self->{dbpass};
}
###
sub dbport {
 my ( $self, $arg ) = @_;
  $self->{dbport} = $arg if (defined $arg);
 return $self->{dbport};
}
###
sub dir {
 my ( $self, $arg ) = @_;
  $self->{dir} = $arg if (defined $arg);
 return $self->{dir};
}
###
sub file {
 my ( $self, $arg ) = @_;
 $self->{file} = $arg if (defined $arg);
 
 return $self->{file};
}
###
sub type {
 my ( $self, $arg ) = @_;
  $self->{type} = $arg if (defined $arg);
 
 return $self->{type};
}
###

sub verbose {
 my ( $self, $arg ) = @_;
  $self->{verbose} = $arg if (defined $arg);
 return $self->{verbose};
}
###
sub descend {
 my ( $self, $arg ) = @_;
  $self->{descend} = $arg if (defined $arg);
 return $self->{descend};
}
###

sub debug {
 my ( $self, $arg ) = @_;
  $self->{debug} = $arg if (defined $arg);
 return $self->{debug};
}
###
sub use_dir {
 my ( $self, $arg ) = @_;
  $self->{use_dir} = $arg if (defined $arg);

 return $self->{use_dir};
}
###

sub list_file {
 my ( $self, $arg ) = @_;
  $self->{list_file} = $arg if (defined $arg);
 return $self->{list_file};
}

1;
