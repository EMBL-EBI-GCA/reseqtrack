package ReseqTrack::Tools::CurrentTreeDiffer;

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::GeneralUtils qw(current_date);
use ReseqTrack::Tools::FileUtils qw(assign_type_by_filename);
use File::Basename qw(fileparse);
use File::Compare qw(compare_text);

=pod

=head1 NAME

ReseqTrack::Tools::CurrentTreeDiffer

=head1 SYNOPSIS

Module for diffing two current trees and creating changelogs in the format used for 1000genomes project.

If other projects require a different output format for a current tree or a changelog
then this module should be used as a base module.  The project-specific child module
should overrun various subroutines.

=head1 Example

  my $tree_differ = ReseqTrack::Tools::CurrentTreeMaker->new(
    -skip_regexes => 'current.tree',
    -dir_to_tree => $dir_to_tree,
    -file_adaptor => $fa,
    -options => \%options,
  );
my $tree_differ = ReseqTrack::Tools::CurrentTreeDiffer->new(
                -old_tree => '/path/to/archive/current.tree',
                -new_tree => '/path/to/staging/current.tree',
                -old_changelog_dir => '/path/to/archive',
                -output_dir   => '/path/to/staging',
                -skip_regexes => ['changelog_details', 'CHANGELOG'],
                -file_type_rules => $ftr_adaptor->fetch_all_in_order,
               );
  if ($tree_differ->quick_diff) {
    $tree_differ->run;
    $tree_differ->print_changelogs;
    $tree_differ->print_errors_to_fh();
  }

=cut

# project-specific child class might override this.
my $NULL = 'NULL'; # used in hashes to mark missing md5

# project-specific child class might add to this list.
sub bad_md5s {
  return {
    'd41d8cd98f00b204e9800998ecf8427e' => 1, # an empty file
    '68b329da9893e34099c7d8ad5cb9c940'=> 1, # a single blank line
  };
}


sub new {
  my ($class, @args) = @_;
  my $self = {};
  bless $self, $class;

  my ( $old_tree, $new_tree, $file_type_rules, $old_changelog_dir, $output_dir, $changelog_name, $skip_regexes, $changelog_details_name)
    = rearrange( [
         qw( OLD_TREE NEW_TREE FILE_TYPE_RULES OLD_CHANGELOG_DIR OUTPUT_DIR CHANGELOG_NAME SKIP_REGEXES CHANGELOG_DETAILS_NAME)
        ], @args);

  $self->old_tree($old_tree);
  $self->new_tree($new_tree);
  $self->file_type_rules($file_type_rules);
  $self->old_changelog_dir($old_changelog_dir);
  $self->output_dir($output_dir);
  $self->changelog_name($changelog_name // 'CHANGELOG');
  $self->changelog_details_name($changelog_details_name // 'changelog_details');
  $self->skip_regexes($skip_regexes);

  #set date in the constructor because cron job runs close to midnight.
  $self->constructor_date(current_date());
  $self->constructor_date_hyphens(current_date(hyphens => 1));

  return $self;
}

=head2 run

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Function  : Calls the subroutines in order which construct the changelog details in memory
              Project-specific child classes can override any of these subroutines.
  Example   : $my_tree_differ->run;

=cut

sub run {
  my ($self) = @_;

  my $old_tree = $self->old_tree;
  throw("no old tree file") if !$old_tree;
  throw("old tree file does not exist $old_tree") if ! -f $old_tree;

  my $new_tree = $self->new_tree;
  throw("no new tree file") if !$new_tree;
  throw("new tree file does not exist $new_tree") if ! -f $new_tree;

  $self->preprocess;
  $self->tree_to_hash($old_tree, 'old');
  $self->tree_to_hash($new_tree, 'new');
  $self->diff_hashes;
  $self->postprocess;
  return;
}

=head2 quick_diff

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Function  : Quickly scans the old and new trees to see if they are significantly different.
              Project-specific child class might override this sub if it uses a different format of the tree file.
  Returntype: Boolean, 1 if old and new files are significantly different.  0 if they are essentially the same.
  Example   : $files_are_different = $my_tree_differ->quick_diff;

=cut

sub quick_diff {
  my ($self) = @_;
  my $old_tree = $self->old_tree;
  throw("no old tree file") if !$old_tree;
  throw("old tree file does not exist $old_tree") if ! -f $old_tree;

  my $new_tree = $self->new_tree;
  throw("no new tree file") if !$new_tree;
  throw("new tree file does not exist $new_tree") if ! -f $new_tree;

  return compare_text($old_tree, $new_tree, sub {
      return 0 if $_[0] eq $_[1];
      return 1 if !$_[0];
      return 0 if (split("\t", $_[0]))[1] eq 'directory';
      return 1;
    } );
}


=head2 preprocess

  Arg [1]   : ReseqTrack::Tools::CurrentTreeDiffer
  Function  : This is called by the run subroutine.  This is an opportunity for
              project-specific child classes to gather any other data from other
              sources before the tree file diffs can be constructed.

=cut

sub preprocess {
  my ($self) = @_;
  return;
}

=head2 postprocess

  Arg [1]   : ReseqTrack::Tools::CurrentTreeDiffer
  Function  : This is called by the run subroutine.  This is an opportunity for
              project-specific child classes to do any cleanup after the tree
              file diffs are constructed

=cut

sub postprocess {
  my ($self) = @_;
  return;
}

=head2 print_details_to_fh

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Arg [2]   : A file handle. Default is *STDOUT
  Arg [3]   : Type, string. This is one of the change types recorded by the diff_hashes
              subroutine e.g. 'withdrawn', 'new', 'moved', 'replacement', or possibly
              something else for a project-specific child class
  Function  : Print the changelog details to this file handle.  Can be called multiple times.
  Example   : $my_tree_differ->print_to_fh($fh, 'replacement')

=cut

sub print_details_to_fh {
  my ($self, $fh, $type) = @_;
  $fh ||= *STDOUT;
  foreach my $change (@{$self->log_change($type)}) {
    print $fh join("\t", @$change), "\n";
  }
  return;
}

=head2 print_details_to_file

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Arg [2]   : A file path
  Arg [3]   : Type, string. This is one of the change types recorded by the diff_hashes
              subroutine e.g. 'withdrawn', 'new', 'moved', 'replacement', or possibly
              something else for a project-specific child class
  Function  : Print the changelog details to this file path.  Can be called multiple times.
  Example   : $my_tree_differ->print_to_fh('/path/to/changelog_details_20140101_replacement', 'replacement')

=cut

sub print_details_to_file {
  my ($self, $output_file, $type) = @_;
  open my $fh, '>', $output_file or throw("could not open $output_file $!");
  $self->print_details_to_fh($fh, $type);
  close $fh;
}

=head2 print_changelogs

  Arg [1]   : ReseqTrack::Tools::CurrentTreeMaker
  Function  : This sub chooses sensible names for the output files and prints them where necessary
              Output files follow the 1000genomes conventions. Project-specific child classes
              should override this subroutine if the changelogs should be handled differently.
  Example   : $my_tree_differ->print_changelogs();

=cut

sub print_changelogs {
  my ($self) = @_;
  my $logged_changes = $self->log_change;
  return if ! keys %$logged_changes;

  my $file_type_rules = $self->file_type_rules;
  throw("this print_changelogs sub requires file_type_rules") if !$file_type_rules;

  my $old_changelog = $self->old_changelog_dir . '/' . $self->changelog_name;
  my $new_changelog = $self->output_dir . '/' . $self->changelog_name;
  my @old_changelog_lines;
  open my $IN, '<', $old_changelog or throw("could not open $old_changelog $!");
  @old_changelog_lines = <$IN>;
  close $IN;

  my $date_hyphens = $self->constructor_date_hyphens();
  my $date = $self->constructor_date();
  my $output_dir = $self->output_dir;
  my @changelog_lines;
  my $changelog_details_name = $self->changelog_details_name;
  push(@changelog_lines, $date_hyphens, "\n", "\n");
  while (my ($change, $details_list) = each %$logged_changes) {
    my %file_types;
    foreach my $details (@$details_list) {
      my $file_type = assign_type_by_filename($details->[0], $file_type_rules);
      $file_types{$file_type} = 1;
    }
    my $short_path = "$changelog_details_name/${changelog_details_name}_${date}_${change}";
    my $full_path = $output_dir . '/' . $short_path;
    $self->print_details_to_file($full_path, $change);
    $self->output_files($full_path);
    push(@changelog_lines, 'Modification to: ' . join(',', map {lc} keys %file_types) . "\n");
    push(@changelog_lines, "\n");
    push(@changelog_lines, "Details can be found in\n");
    push(@changelog_lines, $short_path . "\n\n");
  }

  open my $OUT, '>', $new_changelog or throw("could not open $new_changelog $!");
  print $OUT @changelog_lines;
  print $OUT @old_changelog_lines;
  close $OUT;
  $self->output_files($new_changelog);
}

=head2 print_errors_to_fh

  Arg [1]   : ReseqTrack::Tools::CurrentTreeDiffer
  Arg [2]   : A file handle. Default is *STDERR
  Function  : Print logged error messages to this file handle.  Can be called multiple times.
  Example   : $my_tree_differ->print_errors_to_fh(*STDERR);

=cut

sub print_errors_to_fh {
  my ($self, $fh) = @_;
  $fh ||= *STDERR;
  print $fh map {$_ . "\n"} @{$self->log_error};
}

=head2 tree_to_hash

  Arg [1]   : ReseqTrack::Tools::CurrentTreeDiffer
  Arg [2]   : path to tree file
  Arg [3]   : label, old or new
  Function  : This is called by the run subroutine. Builds hashes of the tree file
              content. These hashes will later get parsed by the diff_hashes
              subroutine.
              Project-specific child classes may want to build hashes
              if the child class also overrides the diff_hashes subroutine
  Example   : $my_tree_differ->tree_to_hash('/path/to/current.tree', 'old')

=cut

sub tree_to_hash {
  my ($self, $tree_file, $label) = @_;
  my $md5_to_path = $self->tree_hashes->{'md5_to_path'};
  my $path_to_md5 = $self->tree_hashes->{'path_to_md5'};
  my $skip_regexes = $self->skip_regexes;
  open my $fh, '<', $tree_file or throw("could not open $tree_file $!");
  LINE:
  while (my $line = <$fh>) {
    chomp $line;
    my ($path, $type, $size, $date, $md5) = split("\t", $line);
    next LINE if $type eq 'directory';
    next LINE if grep {$path =~ /$_/} @$skip_regexes;
    if (!$md5) {
      $self->log_error("no md5 for $path in $tree_file");
      $path_to_md5->{$path}->{$label} = $NULL;
      next LINE;
    }
    push(@{$md5_to_path->{$md5}->{$label}}, $path);
    $path_to_md5->{$path}->{$label} = $md5;
  }
  close $fh;
  $self->tree_hashes('md5_to_path', $md5_to_path);
  $self->tree_hashes('path_to_md5', $path_to_md5);
  return;
}

=head2 diff_hashes

  Arg [1]   : ReseqTrack::Tools::CurrentTreeDiffer
  Function  : This is called by the run subroutine. Hashes built by the tree_to_hash
              subroutine are parsed to identify differences between the new and old hash.
              Project-specific child classes may override this subroutine if the
              changelogs have a different format to 1000genomes
  Example   : $my_tree_differ->diff_hashes();

=cut

sub diff_hashes {
  my ($self) = @_;
  my $md5_to_path = $self->tree_hashes->{'md5_to_path'};
  my $path_to_md5 = $self->tree_hashes->{'path_to_md5'};
  my $bad_md5s = $self->bad_md5s;
  FILE:
  while ( my ($path, $path_hash) = each %$path_to_md5) {
    my ($new_md5, $old_md5) = @{$path_hash}{'new', 'old'};
    if ($new_md5 && $old_md5) {
      if ($new_md5 eq $NULL && $old_md5 eq $NULL) {
        $self->log_error("$path: No md5 in old tree. No md5 in new tree. Ignoring this file");
      }
      elsif ($new_md5 eq $NULL) {
        $self->log_error("$path: Has md5 in old tree. No md5 in new tree. Recording it as replaced file");
        $self->log_change('replacement', $path);
      }
      elsif ($old_md5 eq $NULL) {
        $self->log_error("$path: No md5 in old tree. Has md5 in new tree. Ignoring this file");
      }
      elsif ($new_md5 ne $old_md5) {
        $self->log_change('replacement', $path);
      }
    }
    elsif ($new_md5) {
      my $old_files_match_md5 = $md5_to_path->{$new_md5}->{'old'} // [];
      my $new_files_match_md5 = $md5_to_path->{$new_md5}->{'new'} // [];

      if ( @$old_files_match_md5 == 1 && @$new_files_match_md5 == 1 && $new_md5 ne $NULL && !$bad_md5s->{$new_md5}) {
        if ($path_to_md5->{$old_files_match_md5->[0]}->{'new'}) {
          $self->log_change('new', $path);
        }
        $self->log_change('moved', $old_files_match_md5->[0], $path);
      }
      else {
        my $filename = fileparse($path);
        my @other_new_files_match_filename = grep {fileparse($_) eq $filename}
                              grep {! $path_to_md5->{$_}->{'old'}}
                              grep {$_ ne $path}
                              @$new_files_match_md5;
        my @old_files_match_filename = grep {fileparse($_) eq $filename}
                              grep {! $path_to_md5->{$_}->{'new'}}
                              @$old_files_match_md5;
        if (@other_new_files_match_filename) {
          $self->log_error("$path: Not in old tree. Duplicated file name and md5 in new tree. Recording it as a new file");
          $self->log_change('new', $path);
        }
        elsif (@old_files_match_filename >1) {
          $self->log_error("$path: Not in old tree. Duplicated file name and md5 in old tree. Recording it as a new file");
          $self->log_change('new', $path);
        }
        elsif (@old_files_match_filename == 1) {
          $self->log_change('moved', $old_files_match_md5->[0], $path);
        }
        elsif ($new_md5 eq $NULL) {
          $self->log_error("$path: Not in old tree. No md5 in new tree. Recording it as a new file");
          $self->log_change('new', $path);
        }
        elsif ($bad_md5s->{$new_md5}) {
          $self->log_error("$path: Not in old tree. Bad md5 in new tree. Recording it as a new file");
          $self->log_change('new', $path);
        }
        else {
          $self->log_change('new', $path);
        }
      }
    }
    elsif ($old_md5) {
      my $old_files_match_md5 = $md5_to_path->{$old_md5}->{'old'} // [];
      my $new_files_match_md5 = $md5_to_path->{$old_md5}->{'new'} // [];
      if ( @$old_files_match_md5 == 1 && @$new_files_match_md5 == 1 && $old_md5 ne $NULL && !$bad_md5s->{$old_md5}) {
        if ($path_to_md5->{$new_files_match_md5->[0]}->{'old'}) {
          $self->log_change('withdrawn', $path);
        }
      }
      else {
        my $filename = fileparse($path);
        my @other_old_files_match_filename = grep {fileparse($_) eq $filename}
                              grep {! $path_to_md5->{$_}->{'new'}}
                              grep {$_ ne $path}
                              @$old_files_match_md5;
        my @new_files_match_filename = grep {fileparse($_) eq $filename}
                              grep {! $path_to_md5->{$_}->{'old'}}
                              @$new_files_match_md5;
        if (@other_old_files_match_filename) {
          $self->log_error("$path: Duplicated file name and md5 in old tree. Not in new tree. Recording it as a withdrawn file");
          $self->log_change('withdrawn', $path);
        }
        elsif (@new_files_match_filename >1) {
          $self->log_error("$path: Duplicated file name and md5 in new tree. Not in new tree Recording it as a new file");
          $self->log_change('withdrawn', $path);
        }
        elsif (!@new_files_match_filename) {
          if ($old_md5 eq $NULL) {
            $self->log_error("$path: No md5 in old tree. Not in new tree. Recording it as a new file");
            $self->log_change('withdrawn', $path);
          }
          elsif ($bad_md5s->{$old_md5}) {
            $self->log_error("$path: No md5 in old tree. Bad md5 in new tree. Recording it as a new file");
            $self->log_change('withdrawn', $path);
          }
          else {
            $self->log_change('withdrawn', $path);
          }
        }
      }

    }
  }
  return;
}


sub log_error {
  my ($self, $arg) = @_;
  if ($arg) {
    push(@{$self->{'log_error'}}, $arg);
  }
  return $self->{'log_error'} // [];
}

sub log_change {
  my ($self, $type, @args) = @_;
  if (@_ > 2) {
    push(@{$self->{'log_change'}->{$type}}, [@args]);
  }
  if (@_ == 2) {
    return $self->{'log_change'}->{$type};
  }
  else {
    return $self->{'log_change'} // {};
  }
}

sub old_changelog_dir {
  my ($self, $dir) = @_;
  if (@_ >1) {
    $dir =~ s{//}{/}g;
    $dir =~ s/\/$//;
    $self->{'old_changelog_dir'} = $dir;
  }
  return $self->{'old_changelog_dir'};
}
sub changelog_name {
  my ($self, $arg) = @_;
  if (@_ >1) {
    $self->{'changelog_name'} = $arg;
  }
  return $self->{'changelog_name'};
}
sub changelog_details_name {
  my ($self, $arg) = @_;
  if (@_ >1) {
    $self->{'changelog_details_name'} = $arg;
  }
  return $self->{'changelog_details_name'};
}
sub output_dir {
  my ($self, $dir) = @_;
  if (@_ >1) {
    $dir =~ s{//}{/}g;
    $dir =~ s/\/$//;
    $self->{'output_dir'} = $dir;
  }
  return $self->{'output_dir'};
}
sub output_files {
  my ($self, $file) = @_;
  if ($file) {
    push(@{$self->{'output_files'}}, $file);
  }
  return $self->{'output_files'} // [];
}
sub tree_hashes {
  my ($self, $type, $file) = @_;
  if (@_ > 2) {
    $self->{'tree_hashes'}->{$type} = $file;
  }
  if (@_ == 2) {
    return $self->{'tree_hashes'}->{$type};
  }
  else {
    return $self->{'tree_hashes'} // {};
  }
}

sub file_type_rules {
  my ($self, $ftrs) = @_;
  if ($ftrs) {
    $self->{'file_type_rules'} = $ftrs;
  }
  return $self->{'file_type_rules'};
}

sub old_tree {
  my ($self, $arg) = @_;
  if (@_ >1) {
    $self->{'old_tree'} = $arg;
  }
  return $self->{'old_tree'};
}
sub new_tree {
  my ($self, $arg) = @_;
  if (@_ >1) {
    $self->{'new_tree'} = $arg;
  }
  return $self->{'new_tree'};
}
sub skip_regexes {
  my ($self, $arg) = @_;
  if ($arg) {
    $self->{'skip_regexes'} = ref($arg) eq 'ARRAY' ? $arg : [$arg];
  }
  return $self->{'skip_regexes'} // [];
}

# constructor date is the date when the object was constructed
# (script is run close to midnight)
sub constructor_date {
  my ($self, $arg) = @_;
  if (@_ >1) {
    $self->{'constructor_date'} = $arg;
  }
  return $self->{'constructor_date'};
}
sub constructor_date_hyphens {
  my ($self, $arg) = @_;
  if (@_ >1) {
    $self->{'constructor_date_hyphens'} = $arg;
  }
  return $self->{'constructor_date_hyphens'};
}

1;

