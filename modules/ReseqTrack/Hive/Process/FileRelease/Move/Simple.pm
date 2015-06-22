
package ReseqTrack::Hive::Process::FileRelease::Move::Simple;

use strict;
use File::Basename qw(fileparse);

use base ('ReseqTrack::Hive::Process::FileRelease::Move');

sub derive_path {
  my ($self, $dropbox_path, $file_object) = @_;

  my $derive_path_options = $self->param('derive_path_options');
  my $destination_base_dir = $derive_directory_options->{destination_base_dir};
  throw("this module needs a destination_base_dir") if ! defined $destination_base_dir;

  my ($filename, $dirname) = fileparse($dropbox_path);

  if (my $trim = $derive_directory_options->{trim_dir}) {
    $dirname =~ s/^$trim// or throw("could not trim $dirname $trim");
  }

  my $destination = "$destination_base_dir/$dirname/$filename";
  $destination =~ s/\/\/+/\//g;

  return $destination;

}
1;

