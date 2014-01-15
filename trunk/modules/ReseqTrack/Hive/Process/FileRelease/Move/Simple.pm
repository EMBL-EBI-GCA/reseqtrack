
package ReseqTrack::Hive::Process::FileRelease::Move::Simple;

use strict;
use File::Basename qw(dirname);

use base ('ReseqTrack::Hive::Process::FileRelease::Move');

sub derive_directory {
  my ($self, $dropbox_path, $file_object) = @_;

  my $derive_directory_options = $self->param('derive_directory_options');
  my $destination_base_dir = $derive_directory_options->{destination_base_dir};
  throw("this module needs a destination_base_dir") if ! defined $destination_base_dir;

  my $dirname = dirname($dropbox_path);

  if (my $trim = $derive_directory_options->{trim_dir}) {
    $dirname =~ s/^$trim//;
  }

  my $destination = "$destination_base_dir/$dirname";
  $destination =~ s/\/\/+/\//g;

  return $destination;

}
1;

