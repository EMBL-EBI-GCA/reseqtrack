package ReseqTrack::Tools::Metadata::G1KCheckRunStatus;

use base qw(ReseqTrack::Tools::Metadata::BaseCheckRunStatus);

use strict;
use warnings;

sub skippable_statuses {
  return qw(public);
}

sub collection_types_to_check {
  return qw(FILTERED_FASTQ);
}

sub is_an_ok_file_path {
  my ( $self, $file_path ) = @_;
  return 1 if $file_path =~ /withdrawn/

}

sub correct_path {
  my ( $self, $file_path ) = @_;

  $file_path =~ s/ftp/withdrawn/;

  return $file_path;
}

1;
