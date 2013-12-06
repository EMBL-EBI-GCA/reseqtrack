
package ReseqTrack::Hive::NameFile::BaseNameFile;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw);
use ReseqTrack::Tools::Argument qw(rearrange);
use File::Basename qw(fileparse);

sub method_subs {
  return {
    basic => \&basic_rename,
    default => sub {return},
  };
}


sub new {
  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;

  my ( $method, $params, $files ) = rearrange( [ qw( METHOD PARAMS FILES ) ], @args);

  $self->method($method || 'default');
  $self->params($params);
  $self->files($files);
  
  return $self;
}

sub method {
  my ($self, $method) = @_;
  if ($method) {
    $self->{'method'} = $method;
  }
  return $self->{'method'};
}

sub params {
  my ($self, $params) = @_;
  if ($params) {
    throw "params is not a hashref" if ref($params) ne 'HASH';
    $self->{'params'} = $params;
  }
  return $self->{'params'} // {};
}

sub files {
  my ($self, $files) = @_;
  if ($files) {
    $self->{'files'} = ref($files) eq 'ARRAY' ? $files : [$files];
  }
  return $self->{'files'} // [];
}

sub file_map {
  my ($self, $file_map) = @_;
  if ($file_map) {
    $self->{'file_map'} = $file_map;
  }
  return $self->{'file_map'} // {};
}


sub derive_paths {
  my ($self) = @_;
  my $method = $self->method;
  throw("did not recognise method $method") if ! defined $self->method_subs->{$method};
  return &{ $self->method_subs->{$method}}($self);
}

sub basic_rename {
  my ($self) = @_;
  my $params = $self->params;
  my $files = $self->files;

  my %file_map;

  if (scalar @$files > 1) {
    throw "cannot rename multiple files"
        if $params->{new_name} || $params->{new_basename};
  }

  my $suffixes = $params->{suffix};
  $suffixes = [$suffixes] if !ref($suffixes);

  foreach my $file (@$files) {
    my $new_full_path;
    if ($params->{new_full_path}) {
      $new_full_path = $params->{new_full_path};
    }
    else {
      my ($old_basename, $old_dir, $old_suffix) = fileparse($file, @$suffixes);
      throw "unexpected suffix" if $params->{suffix} && !$old_suffix;
      my $new_dir = $params->{new_dir} // $old_dir;
      my $new_filename;
      if ($params->{new_name}) {
        $new_filename = $params->{new_name};
      }
      else {
        $new_filename = $params->{new_basename} // $old_basename;
        if ($params->{add_datestamp}) {
          my ($year, $month, $day) = (localtime(time))[5,4,3];
          $new_filename .= sprintf(".%04d%02d%02d", $year+1900, $month+1, $day);
        }
        if ($old_suffix) {
          $new_filename .= $old_suffix;
        }
      }
      $new_full_path = "$new_dir/$new_filename";
    }

    $new_full_path =~ s{//+}{/}g;
    $file_map{$file} = $new_full_path;
  }
  $self->file_map(\%file_map);
}

1;
