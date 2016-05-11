package ReseqTrack::Tools::CurrentTreeMaker;

use strict;
use warnings;

use File::Find qw();
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists);
use File::stat;
use Time::localtime;

=pod

=head1 NAME

ReseqTrack::Tools::CurrentTreeMaker

=head1 SYNOPSIS

Module for creating current tree in the format used for 1000genomes project.

If other projects require a different output format in a current tree then this
modules should be used as a base module.  The project-specific child module
should overrun various subroutines.

=head1 Example

  my $tree_maker = ReseqTrack::Tools::CurrentTreeMaker->new(
    -skip_regexes => 'current.tree',
    -dir_to_tree => $dir_to_tree,
    -file_adaptor => $fa,
    -options => \%options,
  );
  $tree_maker->run;
  $tree_maker->print_to_file($new_tree_path);
  $tree_maker->print_errors_to_fh();

=cut

sub new {
  my ($class, @args) = @_;
  my $self = {};
  bless $self, $class;

  my ( $dir_to_tree, $skip_regexes, $trim_dir, $file_adaptor, $files_to_tree, $options)
    = rearrange( [
         qw( DIR_TO_TREE SKIP_REGEXES TRIM_DIR FILE_ADAPTOR FILES_TO_TREE OPTIONS)
        ], @args);

  $self->skip_regexes($skip_regexes);
  $self->trim_dir($trim_dir);
  $self->dir_to_tree($dir_to_tree);
  $self->file_adaptor($file_adaptor);
  $self->files_to_tree($files_to_tree);

  $self->options($self->DEFAULT_OPTIONS);
  $self->options($options) if ($options);

  return $self;
}

=head2 DEFAULT_OPTIONS
  
  child class can override this module to list other project-specific default options

=cut

sub DEFAULT_OPTIONS {return {
    skip_base_directories => 0,# This is here for 1000genomes project only.
    dont_use_nlink => 0
  };}

=head2 run

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Function  : Calls the subroutines in order which construct the current tree in memory
              Project-specific child classes can override any of these subroutines.
  Example   : $my_tree_maker->run;

=cut


sub run {
  my ($self) = @_;
  my $files_to_tree = $self->files_to_tree;
  my $dir_to_tree = $self->dir_to_tree;
  throw("need a files or a dir to tree") if !$dir_to_tree && !@$files_to_tree;
  throw("dir does not exist $dir_to_tree") if $dir_to_tree && ! -d $dir_to_tree;

  $self->make_db_file_hash;
  $self->preprocess;
  if (@$files_to_tree) {
    $self->process_files_in_list;
  }
  if ($dir_to_tree) {
    $self->process_files_in_dir;
  }
  $self->postprocess;
}


=head2 preprocess

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Function  : This is called by the run subroutine.  This is an opportunity for
              project-specific child classes to gather any other data from other
              sources before the tree file can be constructed.

=cut

sub preprocess {
  my ($self) = @_;
  return;
}

=head2 postprocess

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Function  : This is called by the run subroutine.  This is an opportunity for
              project-specific child classes to do any cleanup after the tree
              file is constructed

=cut

sub postprocess {
  my ($self) = @_;
  return;
}

=head2 print_to_fh

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Arg [2]   : A file handle. Default is *STDOUT
  Function  : Print the output to this file handle.  Can be called multiple times.
  Example   : $my_tree_maker->print_to_fh(*STDOUT);

=cut


sub print_to_fh {
  my ($self, $fh) = @_;
  $fh ||= *STDOUT;
  print $fh map {$_ . "\n"} @{$self->tree_lines};
  return;
}

=head2 print_to_file

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Arg [2]   : A file path
  Function  : Print the output to this file path.  Can be called multiple times.
  Example   : $my_tree_maker->print_to_file('/path/to/current.tree');

=cut

sub print_to_file {
  my ($self, $output_file) = @_;
  open my $fh, '>', $output_file or throw("could not open $output_file $!");
  $self->print_to_fh($fh);
  close $fh;
}

=head2 print_errors_to_fh

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Arg [2]   : A file handle. Default is *STDERR
  Function  : Print logged error messages to this file handle.  Can be called multiple times.
  Example   : $my_tree_maker->print_errors_to_fh(*STDERR);

=cut

sub print_errors_to_fh {
  my ($self, $fh) = @_;
  $fh ||= *STDERR;
  print $fh map {$_ . "\n"} @{$self->log_error};
}

=head2 process_files_in_list

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Function  : This is called by the run subroutine.
              Looks at files passed directly to the CurrentTreeMaker object.
              Calls the make_file_line subroutine for each of these files.

=cut

sub process_files_in_list {
  my ($self) = @_;
  my $files_to_tree = $self->files_to_tree;
  my $trim_dir = $self->trim_dir;
  my $skip_regexes = $self->skip_regexes;
  FILE:
  foreach my $full_path (@$files_to_tree) {
    next FILE if grep {$full_path =~ /$_/} @$skip_regexes;
    my $trimmed_path = $full_path;
    if ($trim_dir) {
      $trimmed_path =~ s/^$trim_dir//;
    }
    $self->make_file_line($trimmed_path, $full_path);
  }
  return;
}

=head2 process_files_in_dir

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Function  : This is called by the run subroutine.
              Looks at files and directories in the dir_to_tree
              For each of these it calls the subroutines make_dir_line or make_file_line

=cut

sub process_files_in_dir {
  my ($self) = @_;
  my $dir_to_tree = $self->dir_to_tree;
  my $skip_regexes = $self->skip_regexes;
  my $trim_dir = $self->trim_dir;
  if (!defined $trim_dir) {
    ($trim_dir) = $dir_to_tree =~ /^(.*?)[^\/]+\/*$/;
  }
  my $options = $self->options;
  
  if ($options->{'dont_use_nlink'}){
    $File::Find::dont_use_nlink = 1;
  }

  File::Find::find({
                    wanted => sub{
                        return if -d $_;
                        return if grep {$File::Find::name =~ /$_/} @$skip_regexes;
                        my $trimmed_path = $File::Find::name;
                        $trimmed_path =~ s/^$trim_dir//;
                        $self->make_file_line($trimmed_path, $File::Find::name);
                      },
                    preprocess => sub {
                        return if grep {$File::Find::dir =~ $_} @$skip_regexes;
                        my @dir_content = grep {! /^\./} @_;
                        if ($options->{'skip_base_directories'} && ! grep {-f $_} @dir_content) {
                          return @dir_content;
                        }
                        my $trimmed_path = $File::Find::dir;
                        $trimmed_path =~ s/^$trim_dir//;
                        $self->make_dir_line($trimmed_path, $File::Find::dir);
                        return @dir_content;
                      },
                    }, $dir_to_tree);
  return;
}

=head2 make_file_line

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Arg [2]   : Trimmed file path to output in the the tree file.
  Arg [3]   : Full path of the file on the filesystem
  Function  : This is called once for each file that appears in the output file
              Project-specific child classes should override this subroutine because
              it dictates the format of the final output file

=cut

sub make_file_line {
  my ($self, $trimmed_path, $full_path) = @_;
  my $file_size = -s $full_path;
  my $mtime = ctime(stat($full_path)->mtime);
  my $line = join("\t", $trimmed_path, 'file', $file_size, $mtime);

  my $db_file = $self->db_file_hash->{$full_path};
  if (!$db_file) {
    $self->log_error("nothing in database for $full_path");
    $self->tree_lines($line);
    return;
  }
  my $md5 = $db_file->md5;
  if (!$md5) {
      $self->log_error("no md5 in database for $full_path");
      $self->tree_lines($line);
      return;
  }
  $line .= "\t$md5";
  $self->tree_lines($line);
  return;
}

=head2 make_dir_line

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Arg [2]   : Trimmed dir path to output in the the tree file.
  Arg [3]   : Full dir of the file on the filesystem
  Function  : This is called once for each directory that appears in the output file
              Project-specific child classes should override this subroutine because
              it dictates the format of the final output file

=cut

sub make_dir_line {
  my ($self, $trimmed_path, $full_path) = @_;
  my $dir_size = -s $full_path;
  my $mtime = ctime(stat($full_path)->mtime);
  my $print_path = $full_path;
  $self->tree_lines(join("\t", $trimmed_path, 'directory', $dir_size, $mtime));
  return;
}


=head2 make_db_file_hash

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Function  : This is called by the run subroutine. It builds a hash of files in the database.
              Project-specific child classes should override this subroutine if it
              requires the has to contain any additional information.

=cut

sub make_db_file_hash {
  my ($self) = @_;
  my $fa = $self->file_adaptor or throw("no file adaptor");
  my $dir_to_tree = $self->dir_to_tree;
  my $files_to_tree = $self->files_to_tree;
  my %db_files;
  if ($dir_to_tree) {
    foreach my $file (@{$fa->fetch_by_dirname($dir_to_tree)}) {
      $db_files{$file->name} = $file;
    }
  }
  if (@$files_to_tree) {
    foreach my $path (@$files_to_tree) {
      $db_files{$path} = $fa->fetch_by_name($path);
    }
  }
  $self->db_file_hash(\%db_files);
  return;
}

sub log_error {
  my ($self, $arg) = @_;
  if ($arg) {
    push(@{$self->{'log_error'}}, $arg);
  }
  return $self->{'log_error'} // [];
}

sub tree_lines {
  my ($self, $arg) = @_;
  if ($arg) {
    push(@{$self->{'tree_lines'}}, $arg);
  }
  return $self->{'tree_lines'} // [];
}

sub db_file_hash {
  my ($self, $hash) = @_;
  if ($hash) {
    $self->{'db_file_hash'} = $hash;
  }
  return $self->{'db_file_hash'} // {};
}

sub dir_to_tree {
  my ($self, $dir) = @_;
  if ($dir) {
    $dir =~ s{//}{/}g;
    $dir =~ s/\/$//;
    $self->{'dir_to_tree'} = $dir;
  }
  return $self->{'dir_to_tree'};
}
sub files_to_tree {
  my ($self, $arg) = @_;
  if ($arg) {
    my $files = ref($arg) eq 'ARRAY' ? $arg : [$arg];
    foreach my $file (@$files) {
      $file =~ s{//}{/}g;
    }
    $self->{'files_to_tree'} = $files;
  }
  return $self->{'files_to_tree'} // [];
}
sub trim_dir {
  my ($self, $dir) = @_;
  if (@_ >1) {
    $self->{'trim_dir'} = $dir;
  }
  return $self->{'trim_dir'};
}
sub skip_regexes {
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'skip_regexes'} = ref($arg) eq 'ARRAY' ? $arg : [$arg];
  }
  return $self->{'skip_regexes'} // [];
}
sub file_adaptor {
  my ($self, $fa) = @_;
  if ($fa) {
    $self->{'file_adaptor'} = $fa;
  }
  return $self->{'file_adaptor'};
}

sub file_list {
  my ($self, $list) = @_;
  if ($list) {
    $self->{'file_list'} = $list;
  }
  return $self->{'file_list'};
}

# This options sub looks complicated but it is identical to the one in RunProgram
sub options {
  my ($self, @args) = @_;

  $self->{'options'} ||= {};  
  my $num_of_args = scalar(@args);

  # no arguments, return all the options
  if ($num_of_args == 0){
    return $self->{'options'};
  }
  elsif ($num_of_args == 1){
    my $ref_type = ref($args[0]);

    if (! $ref_type){
      # arg is a scalar
      return $self->{'options'}->{$args[0]};
    }
    elsif ($ref_type eq 'HASH'){
      # merge these options with any existing

      while (my ($name, $value) = each %{$args[0]}) {
        $self->{'options'}->{$name} = $value;
      }

      return $self->{'options'};
    }
    else {
      throw("Cannot set options with a $ref_type reference");
    }

  }
  else{
    my ($option_name,$option_value) = @args;
    throw( "option_name not specified") if (! $option_name);
    $self->{'options'}->{$option_name} = $option_value;
    return $self->{'options'}->{$option_name};
  }
}

1;

