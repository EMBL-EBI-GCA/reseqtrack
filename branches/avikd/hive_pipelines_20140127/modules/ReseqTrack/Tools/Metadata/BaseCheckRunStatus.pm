package ReseqTrack::Tools::Metadata::BaseCheckRunStatus;
use strict;
use warnings;
use base qw(ReseqTrack::Tools::Metadata::BaseMetadataAddIn);

use ReseqTrack::Tools::Exception qw(throw);

=pod

=head1 NAME

ReseqTrack::Tools::Metadata::BaseCheckRunStatus

=head1 SYNOPSIS

This module checks that files associated with withdrawn runs are in the correct location.
It should be extended per project - see G1kCheckRunStatus for an example of this.

Runs will be checked unless they have a whitelisted status. For G1K, this should be public. Other projects may need to add private.

Collections will be found using the run_source_id (e.g. ERR123456) and the types listed in $self->collection_types_to_check

Each file 

These methods must be implemented by classes extending this one:

skippable_statuses
collection_types_to_check
is_an_ok_file_path
correct_path

=head1 OPTIONS

=cut

sub check_run {
  my ( $self, $run, $current_copy ) = @_;
  my $collection_adaptor = $self->reseq_db->get_CollectionAdaptor();

  if ( !$self->is_a_skippable_status( $run->status ) ) {
    for my $collection_type ( $self->collection_types_to_check() ) {
      my $collection =
        $collection_adaptor->fetch_by_name_and_type( $run->run_source_id,
        $collection_type );
      next unless $collection;
      throw(
"Found non-file collection for $run->run_source_id $collection_type, cannot check paths"
      ) if ( $collection->table_name ne 'file' );

      for my $file ( @{$collection->others} ) {
        if ( !$self->is_an_ok_file_path( $file->name ) ) {
          my $new_path = $self->correct_path( $file->name );
          $self->file_hash( $file->name, $new_path );
        }
      }
    }

  }
}

sub report {
  my ($self) = @_;

  my $file_hash    = $self->file_hash();
  my $status_count = keys(%$file_hash);
  my $log_fh       = $self->log_fh();

# Delierate use of STDOUT instead of $log_fh - the log_fh should contain the file list, warning about the number of errors should go to stdout
  print STDOUT "There are $status_count status issues to resolve$/"
    if ($status_count);

  while ( my ( $original_path, $suggested_path ) = each %$file_hash ) {
    print $log_fh "$original_path\t$suggested_path$/";
  }
}

=head2 correct_path

  Arg [file_path]   :
      path of file that is in the wrong location based on the run status
  
  Function  : Creates a new ReseqTrack::Tools::RunProgram object 
  Returntype: String, suggested new location for this file
	Implementation: See G1kCheckRunStatus for an example implementation

=cut

sub correct_path {
  my ( $self, $file_path ) = @_;
  throw("implement me!");
}

sub is_a_skippable_status {
  my ( $self, $status ) = @_;

  my @match = grep { $_ eq $status } $self->skippable_statuses();

  return (@match) ? 1 : undef;
}

=head2 is_an_ok_file_path
  
  Function  : Is this file path correct for a withdrawn file?
  Returntype: 1 for true, undef for false
	Implementation: See G1kCheckRunStatus for an example implementation

=cut

sub is_an_ok_file_path {
  my ( $self, $file_path ) = @_;
  throw("implement me!");
}

=head2 skippable_statuses
  
  Function  : Give a list of run statuses that can be ignored
  Returntype: List of strings
	Implementation: See G1kCheckRunStatus for an example implementation

=cut

sub skippable_statuses {
  throw("implement me!");
}

=head2 collection_types_to_check

  Arg [file_path]   :
      path of file that is in the wrong location based on the run status
  
  Function  : Creates a new ReseqTrack::Tools::RunProgram object 
  Returntype: String, suggested new location for this file
	Implementation: See G1kCheckRunStatus for an example implementation

=cut

sub collection_types_to_check {
  throw("implement me!");
}

sub file_hash {
  my ( $self, $k, $v ) = @_;

  if ( $k && $v ) {
    $self->{file_hash}->{$k} = $v;
  }

  return $self->{file_hash};
}

1;
