
=head1 NAME
ReseqTrack::Tools::Loader

=head1 SYNOPSIS

Misguided attempt to get Loader stuff to work
=head1 Example


=cut

package ReseqTrack::Tools::Loader;
use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileUtils
  qw(create_objects_from_path_list create_history assign_type check_type);
use ReseqTrack::Tools::HostUtils qw(get_host_object);
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
      $verbose, $descend,   $debug
   )
   = rearrange(
  [
   qw(
     DB     DBHOST     DBNAME     DBUSER     DBPASS
     DBPORT     DIR     FILE     TYPE  LIST_FILE
     VERBOSE     DESCEND     DEBUG )
  ],
  @args
   );




  #defaults
 $self->type("MUST_FIX")  ;
 $self->descend('1') ;

 
 

 $self->dbhost($dbhost);
 $self->dbname($dbname);
 $self->dbuser($dbuser);
 $self->dbpass($dbpass);
 $self->dbport($dbport);
# $self->hostname($hostname);

 if ($db) {
   $self->db($db);
 }

 $self->dir($dir);
 $self->file($file);
 $self->type($type);
 $self->list_file($list_file);
 $self->verbose($verbose);
 $self->descend($descend);
 $self->debug($debug);
  return $self;
}
###

sub create_objects{
  my ($self) = @_;
  throw("Must implement a 'create_objects method' in ".$self);
}

sub load_objects{
  my ($self) = @_;
  throw("Must implement a 'load_objects method' in ".$self);
}

sub sanity_check_objects{
   my ($self) = @_;
   throw("Must implement a 'sanity_check_objects' in ".$self);
}

sub process_input{
   my ($self) = @_;
   throw("Must implement a 'process_input' in ".$self);
}

####
sub create_DBAdaptor {
 my $self = shift;
 my $db =
   ReseqTrack::DBSQL::DBAdaptor->new(
				     -host   => $self->dbhost,
				     -user   => $self->dbuser,
				     -port   => $self->dbport,
				     -dbname => $self->dbname,
				     -pass   => $self->dbpass,
				    );
 $db->dbc->disconnect_when_inactive(1);
 throw("Failed to create a database adaptor using ".$self->dbhost." ".$self->dbname)
   unless($db);
 return $db;

}

sub add_files_from_dir {
 my $self         = shift;
 my $descend = $self->descend;
 my @paths;

 return if ( !$self->dir );
 
 if ($self->dir) {
  my ( $paths, $hash ) = list_files_in_dir( $self->dir, 1 );
  unless ($descend) {
    my $dir = $self->dir;
    $dir =~ s/\/$//s;
    $paths = $hash->{$self->dir};
  }
  throw( "No paths were found in " . $self->dir ) if ( !$paths || @$paths == 0 );
  push( @paths, @$paths ) if ( $paths && @$paths >= 1 );
 }

 $self->file_paths( \@paths ) ;
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
 
 $self->file_paths( \@paths ) ;
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
 $self->file_paths( \@paths ) ;

}
####

 
 sub file_paths{
  my ($self, $arg) = @_;
 
  
  if($arg){
    if(ref($arg) eq 'ARRAY'){
      push(@{$self->{file_paths}}, @$arg);
    }else{
      push(@{$self->{file_paths}}, $arg);
    }
  }
  return $self->{file_paths}
}
 
###
sub db {
 my ( $self, $arg ) = @_;

  $self->{db} = $arg if ( $arg) ;
  
   if(!$self->{db}){
    # print "Have to create DBAdaptor\n";
    $self->{db} = $self->create_DBAdaptor();
  }
 return $self->{db};
}
###
sub dbhost {
 my ( $self, $arg ) = @_;
   $self->{dbhost} = $arg if ( $arg);
 return $self->{dbhost};
}
###
sub dbname {
 my ( $self, $arg ) = @_;
   $self->{dbname} = $arg if ( $arg);
  return $self->{dbname};
}
###
sub dbuser {
 my ( $self, $arg ) = @_;
   $self->{dbuser} = $arg if ( $arg);
  return $self->{dbuser};
}
###
sub dbpass {
 my ( $self, $arg ) = @_;
  $self->{dbpass} = $arg if ( $arg);
 return $self->{dbpass};
}
###
sub dbport {
 my ( $self, $arg ) = @_;
  $self->{dbport} = $arg if ($arg);
 return $self->{dbport};
}
###

#May move from Loader/File.pm to here. RES. 
#sub hostname {
# my ( $self, $arg ) = @_;
#  $self->{hostname} = $arg if ($arg);
# return $self->{hostname};
#}
###
sub dir {
 my ( $self, $arg ) = @_;
  $self->{dir} = $arg if ( $arg);
 return $self->{dir};
}
###
sub file {
 my ( $self, $arg ) = @_;
 $self->{file} = $arg if ( $arg);
 
 return $self->{file};
}
###
sub type {
 my ( $self, $arg ) = @_;
  $self->{type} = $arg if ( $arg);
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
  $self->{use_dir} = $arg if ( $arg);
 return $self->{use_dir};
}
###

sub list_file {
 my ( $self, $arg ) = @_;
  $self->{list_file} = $arg if ( $arg);
 return $self->{list_file};
}


1;
