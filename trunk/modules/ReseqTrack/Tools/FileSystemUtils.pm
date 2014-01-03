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
use File::Path qw(remove_tree make_path);
use Digest::MD5::File qw(file_md5_hex);
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
  create_tmp_process_dir  );

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
 throw( "Must pass list_files_in_dir a directory not " . $dir ) if ( !-d $dir );

 my %hash;
 my @files;
 File::Find::find( sub{
      return if /^\./;
      return if -d $_;
      push(@files, $File::Find::name);
      my $filename = $get_full_paths ? $File::Find::name : $_;
      push(@{$hash{$File::Find::dir}}, $filename);
  }, $dir);

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
 throw( $file . " must exist in order to be read" ) if ! -e $file;
 my @return;
 open( my $fh, '<', $file ) || throw("Cannot open the file $file $!");

 while (my $line = <$fh>) {
  chomp($line);
  next if ( $line =~ /^\s*$/ );
  push @return, $line;
 }
 close($fh) or throw( "Failed to close $file $!");
 warning( $file . " contained no lines" ) if ( @return == 0 );
 return \@return;
}

=head2 get_md5hash

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
 throw( $file . " must exist in order to be read" ) if ! -e $file;
 open( my $fh, '<', $file ) or throw( "Failed to open $file $!");
 my %filehash;
 my %md5hash;
 while (<$fh>) {
  chomp;
  my ( $md5, $filename ) = split;
  throw("line does not contain md5 and file") if !$filename || ! length($filename);
  throw("md5 is not valid: $md5") if $md5 !~ /^[a-f0-9]{32}$/;

  my $name = $filename;
  $name = basename($filename) if ($trim);
  warning( $name . " already existing in hash" )
    if ( $filehash{$name} );

  $filehash{$name} = $md5;
  $md5hash{$md5}   = $filename;
 }
 close($fh) or throw( "Failed to close $file $!" );
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
 open( my $fh, '<', $file ) or throw( "Failed to open " . $file . " $!" );
 my %hash;
 while (my $string = <$fh>) {
  chomp $string;
  my @values = split("\t", $string);
  push( @{ $hash{ $values[$column] } }, $string );
 }
 close $fh;
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
 open( my $fh, '<', $file ) or throw( "Failed to open " . $file . " $!" );
 my %hash;
 while (my $string = <$fh>) {
  chomp $string;
  my @values = split("\t", $string);
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
 close $fh;
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
 my ( $file ) = @_;
 throw("Can't run md5sum on a non existent file $file") if ( !-e $file );
 my $md5 = file_md5_hex($file);
 return $md5;
}

=head2 find_file

  Arg [1]   : string, pattern/regex
  Arg [2]   : string, directory to start search in
  Arg [3]   : binary flag verbosity
  Function  : run unix find
  Returntype: arrayref of paths which match pattern
  Exceptions: throws if can't run find command
  Example   : my $fastqs = find('.*\.fastq.gz', "/nfs/1000g-work/");

=cut

sub find_file {
 my ( $name, $dir, $verbose ) = @_;
 my @paths;
 File::Find::find( sub{ push(@paths, $File::Find::name) if ( /$name/ ); }, $dir);
 return \@paths;
}

=head2 check_md5

  Arg [1]   : filepath to file
  Arg [2]   : md5 to check against
  Function  : runs md5sum to see if it equals the given md5
  actual files
  Returntype: 0/1 
  Example   : 

=cut

sub check_md5 {
 my ( $file, $md5 ) = @_;
 my $calculated_md5 = run_md5($file);
 return 1 if $md5 eq $calculated_md5;
 return 0;
}

sub dump_dirtree_summary{
  my ($input_dir, $output_file, $skip_regex, $fa, $file_list, $output_prefix) = @_;
  my $no_md5s = 0;
  $skip_regex ||= 'current.tree';

  my ($files, $hash);# = list_files_in_dir($input_dir, 1);

  #provide file list from 'find $input_dur > file.list'
  if ($file_list){
    $files = get_lines_from_file($file_list);
  }
  else{
    ($files) = list_files_in_dir($input_dir, 1);
  }


  my $fh;
  if($output_file){
    open($fh, ">".$output_file) or throw("Failed to open ".$output_file." $!");
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
  $trim =~ s/\/\/+/\//g;
  foreach my $file(@$files){
    next if($file =~ /$skip_regex/);
    my $dir = dirname($file);
    my $label;
    my $md5;
    if($fa){
      $md5 = $file_md5s{$files};
    }
    my $mod_dir = $dir;
    $mod_dir =~ s/$trim\/*//;
    $mod_dir = $mod_dir ? $output_prefix . '/' . $mod_dir : $output_prefix;
    $mod_dir =~ s{//+}{/}g;
    unless($dirs{$mod_dir}){
      my $dir_size = -s $dir;
      my $dir_stamp = ctime(stat($dir)->mtime);
      $label = 'directory';
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
    $file = $output_prefix . '/' . $file;
    $file =~ s{//+}{/}g;
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

  remove_tree($dir, {error => \my $err, verbose=>$verbose});
  if (@$err) {
    my @error_strings;
    foreach my $err_hash (@$err) {
      my ($err_file, $message) = %$err_hash;
      push(@error_strings, "Error creating $err_file: $message.");
    }
    throw join(' ', @error_strings);
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
      	print "$file does not exist. Skipping delete\n" if $verbose;
      	return 1;
      }
    }

    my $result = unlink($file);

    # no need to throw if file has been delete by another process
    if (!$result && -e $file) {
      throw( "file not deleted: $file $!");
    }

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
  throw("$dir is not a directory") if (-e $dir);

  make_path($dir, {mode => 0775, error => \my $err});
  if (@$err) {
    my @error_strings;
    foreach my $err_hash (@$err) {
      my ($err_dir, $message) = %$err_hash;
      push(@error_strings, "Error creating $err_dir: $message.");
    }
    throw join(' ', @error_strings);
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
  if ( $file !~ /^\// ){
    warning("Not full path to $file. But it exists");
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
 
  if ( $file !~ /^\// ){
    warning("Not full path to $file. But it does not exist");
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
  throw "executable not given" if !$executable;

  if ($executable =~ m{/}) {
    throw "executable does not exist: $executable" if (! -e $executable);
    throw "executable is not executable: $executable" if (! -x $executable);
  }
  elsif (! first {-x $_} map {$_.'/'.$executable} @PATH) {
      throw "cannot find executable in path: $executable";
  }
  return 1;
}

=head2 create_tmp_process_dir

  Arg [1]   : path to parent directory in which to make temp dir
  Arg [2]   : basename of the tempdirectory (can be undefined)
  Arg [3]   : cleanup 0/1.  Directory can  be automatically deleted when it goes out of scope in perl
  Function  : checks that the file is an executable
  Returntype: path to temp directory
  Exceptions: throws if no parent directory is given
              throws if can't create a temp directory
              throws if can't chmod the temp directory
  Example   : create_tmp_process_dir('/path/to/dir', 'temp', 1);

=cut


sub create_tmp_process_dir {
  my ($parent_dir, $basename, $cleanup) = @_;
 
  throw "No parent directory specified" if ( ! $parent_dir );
  my $template = $basename ? "$basename." : '';
  $template .= 'XXXXXXXXXX';
  
  my $temp_dir;
  eval{$temp_dir = tempdir( $template, DIR => $parent_dir , CLEANUP => $cleanup);};
  if (!$temp_dir) {
      throw("error creating a temp directory in $parent_dir: $!");
  }

  print "created temp dir $temp_dir \n";

  chmod 0775, $temp_dir or throw("could not chmod $temp_dir $!");
  return $temp_dir;
}

1;
