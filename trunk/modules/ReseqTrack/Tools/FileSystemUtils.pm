package ReseqTrack::Tools::FileSystemUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::File;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use File::Copy;
use File::Basename;
use File::Find ();
use File::stat;
use File::Path;
use File::Temp qw/ tempfile tempdir /;
use Time::localtime;
use List::Util qw (first);
use Env qw( @PATH );

use vars qw (@ISA  @EXPORT);

@ISA    = qw(Exporter);
@EXPORT = qw(
  list_files_in_dir
  get_lines_from_file
  get_md5hash
  get_hash_from_tab_file
  get_columns_from_tab_file
  run_md5
  find_file
  check_md5
  dump_dirtree_summary
  delete_directory
  delete_file
  check_file_exists
  check_file_does_not_exist
  check_directory_exists
  check_executable
  make_directory
  create_tmp_process_dir  );

=head2 get_filenames

  Arg [1]   : string, directory path
  Function  : returns a list of files in given directory, this is a utility method
  to descend directory trees for list_files_in_dir
  Returntype: arrayref of filepaths
  Exceptions: none
  Example   : my $list = get_filenames($dir);

=cut

sub get_filenames {
 my ($dir) = @_;
 no strict;
 sub find(&@) { &File::Find::find }
 *name = *File::Find::name;
 my @files;
 find { push( @files, $name ) } $dir;
 use strict;
 return \@files;
}

=head2 list_files_in_dir

  Arg [1]   : String, directory path
  Function  : get a list of the filenames in the path including full
  descent of the directory tree
  Returntype: arrayref of filenames and hash keys on dir, each containing list
  of files
  Exceptions: throws if not passed a directory, if directory 
  contains ~ files or if it can't open or close the dir
  Example   : my ($files) = list_files_in_dir($dir);

=cut

sub list_files_in_dir {
 my ( $dir, $get_full_paths ) = @_;
 my $name = $File::Find::name;
 throw( "Must pass list_files_in_dir a directory not " . $dir )
   unless ( -d $dir );
 my $list = get_filenames($dir);

 #find { push(@list, $name)} $dir;
 my %hash;
 my @files;
 foreach my $element (@$list) {
  next if ( !$element );
  if ( -d $element ) {
   if ( !$hash{$element} ) {
    $hash{$element} = [];
   }
  }
  else {
   my $file;
   my $dir_path = dirname($element);
   if ($get_full_paths) {
    $file = $element;
   }
   else {
    $file = basename($element);
   }
   my $name = basename($element);
   next if ( $name =~ /^\./ );
   push( @files, $file );
   if ( !$hash{$dir_path} ) {
    $hash{$dir_path} = [];
   }
   push( @{ $hash{$dir_path} }, $file );
  }
 }

 return undef if ( @files == 0 );
 return ( \@files, \%hash );
}

=head2 get_lines_from_file

  Arg [1]   : String, name of file
  Function  : open file and read lines
  Returntype: arrayref of lines from file
  Exceptions: throws if file doesn't exist or for failures 
  in open and closes
  Example   : my @lines = @{get_lines_from_file($file)};

=cut

sub get_lines_from_file {
 my ($file) = @_;
 unless ( -e $file ) {
  throw( $file . " must exist in order to be read" );
 }
 my @return;
 open( FILE, $file ) || throw("Cannot open the file $file");

 foreach my $line (<FILE>) {
  chomp($line);
  next if ( $line =~ /^\s+$/ );
  push @return, $line;
 }
 close(FILE) or throw( "Failed to close " . $file );
 warning( $file . " contained no lines" ) if ( @return == 0 );
 return \@return;
}

=head2 get_Md5hash

  Arg [1]   : String, filename
  Function  : parse file to produce hash of filenames
  pointing to md5 values. Files need to be in format
  md5 filename
  Returntype: hashref, filename->md5
  Exceptions: throws if not given a file, if fails to
  open or close file and if the md5 is zero or the md5 or
  filename is repeated
  Example   : my $md5_hash = get_Md5hash($md5_file);

=cut

sub get_md5hash {
 my ( $file, $trim ) = @_;
 unless ( -e $file ) {
  throw( $file . " must exist in order to be read" );
 }
 open( FH, $file ) or throw( "Failed to open " . $file );
 my %filehash;
 my %md5hash;
 while (<FH>) {
  chomp;
  my ( $md5, $filename ) = split;
  throw( $filename . " has zero md5 in " . $file )
    if ( !$md5 || ( $md5 =~ /^\d+$/ && $md5 == 0 ) );
  warning( $filename . " already existing in hash" )
    if ( $filehash{$filename} );

  #warning($md5." is a duplicate") if($md5hash{$md5});
  if ( !$md5 || length($md5) == 0 ) {
   print $_. " doesn't have an md5\n";
  }
  if ( !$filename || length($filename) == 0 ) {
   print $_. " doesn't have a filename\n";
  }
  my $name = $filename;
  $name = basename($filename) if ($trim);
  $filehash{$name} = $md5;
  $md5hash{$md5}   = $filename;
 }
 close(FH) or throw( "Failed to close " . $file );
 return \%filehash;
}

=head2 get_hash_from_tab_file

  Arg [1]   : string, filepath
  Arg [2]   : int, column index (starting at 0)
  Function  : return a hash keyed on the specified column. The key is the column
  value, the value is an arrayref of the lines associated with that value
  Returntype: hashref
  Exceptions: none
  Example   : my $file_hash = get_hash_from_tab_file($file, 0);

=cut

sub get_hash_from_tab_file {
 my ( $file, $column ) = @_;
 open( FH, $file ) or throw( "Failed to open " . $file . " $!" );
 my %hash;
 while (<FH>) {
  chomp;
  my $string = $_;
  my @values = split /\t/, $string;
  if ( !$hash{ $values[$column] } ) {
   $hash{ $values[$column] } = [];
  }
  push( @{ $hash{ $values[$column] } }, $string );
 }
 return \%hash;
}

=head2 get_columns_from_tab_file

  Arg [1]   : string, file path
  Arg [2]   : int, column of tab delimited string to use a key
  Arg [3]   : int, column of tab delimited string to use a value
  Function  : returns hash based on two columns of a tab file
  Returntype: hashref
  Exceptions: warns if a key already has a value but it will be replaced
  Example   : my $md5_hash = get_columns_from_tab_file($file, 0, 1);

=cut

sub get_columns_from_tab_file {
 my ( $file, $key_column, $value_column ) = @_;
 open( FH, $file ) or throw( "Failed to open " . $file . " $!" );
 my %hash;
 while (<FH>) {
  chomp;
  my $string = $_;
  my @values = split /\t/, $string;
  if ( $hash{ $values[$key_column] } ) {
   warning(   "ReseqTrack::Tools::FileUtils::get_columns_from_tab_file "
            . "assumes that the contents of "
            . $key_column . " in "
            . $file . " are "
            . "unique but they appear not to be, this will means "
            . $value_column
            . "entries in the values of the hash maybe misleading" );

  }
  $hash{ $values[$key_column] } = $values[$value_column];
 }
 return \%hash;
}

=head2 run_md5

  Arg [1]   : string, filepath
  Arg [2]   : string, program path
  Function  : runs md5sum on file
  Returntype: string, 32bit md5 checksum
  Exceptions: throws if file doesn't exist or if commandline fails
  Example   : my $md5 = run_md5($file);

=cut

sub run_md5 {
 my ( $file, $program ) = @_;
 $program = 'md5sum' if ( !$program );
 throw("Can't run md5sum on a non existent file $file") unless ( -e $file );
 my $cmd = $program . " " . $file;
 open( FH, $cmd . " | " ) or throw( "Failed to open " . $cmd );
 my $md5;
 while (<FH>) {
  chomp;
  my @values = split;
  $md5 = $values[0];
 }
 close(FH);
 return $md5;
}

=head2 find_file

  Arg [1]   : string, name/pattern
  Arg [2]   : string, directory to start search in
  Arg [3]   : binary flag verbosity
  Function  : run unix find
  Returntype: arrayref of paths which match pattern
  Exceptions: throws if can't run find command
  Example   : my $fastqs = find("*fastq.gz", "/nfs/1000g-work/");

=cut

sub find_file {
 my ( $name, $dir, $verbose ) = @_;
 my $cmd = "find $dir -name \"" . $name . "\" -print";
 print STDERR $cmd . "\n" if ($verbose);
 open( CMD, $cmd . " | " ) or throw( "Failed to open " . $cmd . " $!" );
 my @paths;
 while (<CMD>) {
  chomp;
  push( @paths, $_ );
 }
 close(CMD);
 return \@paths;
}

=head2 check_md5

  Arg [1]   : filepath to file containing md5
  Arg [2]   : path to md5 program, defaults to md5sum
  Function  : runs md5sum -c with given file and checks md5s in file match md5s of
  actual files
  Returntype: 0/1 
  Exceptions: throws if failed to run cmd
  Example   : 

=cut

sub check_md5 {
 my ( $file, $hash ) = @_;
 my $name           = basename($file);
 my $md5            = $hash->{$name};
 my $calculated_md5 = run_md5($file);
 unless ( $md5 eq $calculated_md5 ) {
  return 0;
 }
 else {
  return 1;
 }
 print STDERR "Shouldn't of reached this point in FileUtils check_md5\n";
}

sub dump_dirtree_summary{
  my ($input_dir, $output_file, $skip_regex, $fa, $file_list) = @_;
  my $no_md5s = 0;
  $skip_regex ||= 'current.tree';

  my ($files, $hash);# = list_files_in_dir($input_dir, 1);

  #provide file list from 'find $input_dur > file.list'
  if ($file_list){
    $files = get_lines_from_file($file_list);
  }
  else{
    ($files, $hash) = list_files_in_dir($input_dir, 1);
  }


  my $fh;
  if($output_file){
    open(FH, ">".$output_file) or throw("Failed to open ".$output_file." $!");
    $fh = \*FH;
  }else{
    $fh = \*STDOUT;
  }
  my %dirs;
  my %file_md5s;
  if($fa){
    my $file_objects = $fa->fetch_all_like_path($input_dir);
    foreach my $file_object(@$file_objects){
      my $md5 = $file_object->md5;
      $md5 = "................................" unless($md5);
      $file_md5s{$file_object->name} = $md5;
    }
  }
  my $trim = $input_dir;
  $trim =~ s/ftp//;
  $trim =~ s/\/\/$/\//;
  foreach my $file(@$files){
    next if($file =~ /$skip_regex/);
    my $dir = dirname($file);
    my $label;
    my $md5;
    if($fa){
      $md5 = $file_md5s{$files};
    }
    my $mod_dir = $dir;
    $mod_dir =~ s/$trim//;
    unless($dirs{$mod_dir}){
      my $dir_size = -s $dir;
      my $dir_stamp = ctime(stat($dir)->mtime);
      $label = 'directory';
      $dir =~ s/$input_dir//;   
      print $fh join("\t", $mod_dir, $label, $dir_size, $dir_stamp);
      print $fh "\t " if($fa);
      print $fh "\n";
      $dirs{$mod_dir} = 1;
    }

    my $md5sum = '';
    $md5sum = $file_md5s{$file};

    #warning($file." has no md5") unless($md5sum);
    if (!$md5sum){
      print STDERR  "$file has no md5\n";
      $no_md5s++;
    }

    my $size = -s $file;
    my $date_string = ctime(stat($file)->mtime);
    $label = 'file';
    $file =~ s/$trim//;
    print $fh join("\t", $file, $label, $size, $date_string);

   # print $fh "\t".$md5sum if($fa);
   # print $fh "\n";

    if  ($md5sum){
      print $fh "\t". $md5sum;
    }
    else{
      print $fh "\t";
    }
    print $fh "\n";
  }
  close($fh);
}

=head2 delete_directory

  Arg [1]   : path to directory.  This directory must exist.
  Function  : deletes the directory and all sub-directories and files.
  Returntype: 0/1 
  Exceptions: throws if deletion fails
  Example   : delete_directory('/path/to/directory');

=cut

sub delete_directory {
  my ($dir, $verbose) = @_;

  eval { rmtree($dir, $verbose)};
  if ($@) {
      throw ($@);
  }
  return 1;
}

=head2 delete_file

  Arg [1]   : path to file.  This file must exist.
  Function  : deletes the file
  Returntype: 0/1 
  Exceptions: throws if deletion fails
  Example   : delete_file('/path/to/file');

=cut

sub delete_file {
    my ($file, $verbose) = @_;

    throw ("no file name") if (!$file);
    print "Deleting $file\n" if ($verbose);

    if ( ! -e $file ){
      if ( ! -l $file ) {
      	print "$file does not exist. Skipping delete\n";
      	return 1;
      }
    }

    unlink($file)
        or throw( "file not deleted: $file $!");

    return 1;
}

=head2 check_directory_exists

  Arg [1]   : path to directory
  Function  : checks that the directory exists
  Returntype: 0/1 
  Exceptions: throws if directory cannot be created
  Example   : check_directory_exists('/path/to/directory');

=cut

sub check_directory_exists {
  my $dir = shift;

  return 1 if (-d $dir);

  eval { mkpath($dir, 0, 0775)};
  if ($@) {
      throw ("Failed to create $dir: $@");
  }
  return 1;
}

=head2 check_file_exists

  Arg [1]   : path to file
  Function  : checks that the file exists and is not a directory
  Returntype: 0/1 
  Exceptions: throws if file does not exist or if it is a directory.
  Example   : check_file_exists('/path/to/file');

=cut

sub check_file_exists {
  my $file = shift;
 
 throw( "$file is a directory, not a file") if ( -d $file );
 throw( "$file does not exist")     if ( !-e $file );
 if ( !( $file =~ /^\// ) ){
  warn "Not full path to $file. But it exists";
  }
  return 1;
}

=head2 check_file_does_not_exist

  Arg [1]   : path to file
  Function  : checks that the file does not exist
  Returntype: 0/1 
  Exceptions: throws if file exists.
  Example   : check_file_does_not_exist('/path/to/file');

=cut

sub check_file_does_not_exist {
  my $file = shift;

  throw("$file already exists") if (-e $file);
 
  if ( !( $file =~ /^\// ) ){
    warn "Not full path to $file. But it does not exist";
  }
  return 1;
}

=head2 check_executable

  Arg [1]   : path to file, or name of executable in path
  Function  : checks that the file is an executable
  Returntype: 0/1 
  Exceptions: throws if file is not an executable
  Example   : check_file_exists('/path/to/executable');
  Example   : check_file_exists('executable_in_path');

=cut

sub check_executable {
  my $executable = shift;

  if ($executable =~ m{/}) {
    throw "executable does not exist: $executable" if (! -e $executable);
    throw "executable is not executable: $executable" if (! -x $executable);
  }
  elsif (! first {-x $_} map {$_.'/'.$executable} @PATH) {
      throw "cannot find executable in path: $executable";
  }
  return 1;
}

=head2 make_directory

  Arg [1]   : path to new directory
  Function  : creates new directory if it does not already exist
  Returntype: 0/1 
  Exceptions: throws if failed to create directory
  Example   : make_directory('/path/to/directory');

=cut

sub make_directory {
    my $dir = shift;

    if ( ! -d $dir) {
        eval { mkpath($dir, 0, 0775)};
        if ($@) {
            throw ($@);
        }
    }

    return 1;
}


sub create_tmp_process_dir {
  my ($parent_dir) = @_;
 
  throw "No parent directory specified" if ( ! $parent_dir );
  
  my $temp_dir = tempdir( DIR => $parent_dir );

  print "processing in would be $temp_dir \n";

  `chmod -R  775 $temp_dir`;
  return ($temp_dir);
}

1;
