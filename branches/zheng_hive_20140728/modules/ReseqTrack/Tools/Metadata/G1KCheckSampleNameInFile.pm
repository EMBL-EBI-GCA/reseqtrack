package ReseqTrack::Tools::Metadata::G1KCheckSampleNameInFile;

use strict;
use warnings;
use base qw(ReseqTrack::Tools::Metadata::BaseMetadataAddIn);

use ReseqTrack::Tools::Exception qw(throw);

=pod

=head1 NAME

ReseqTrack::Tools::Metadata::G1CheckRunStatus

=head1 SYNOPSIS

This module checks that files associated with a run have the sample_alias in the path
Collections will be found using the run_source_id (e.g. ERR123456) and the types listed in $self->collection_types_to_check

=cut

sub check_run {
  my ( $self, $run, $current_copy ) = @_;

  my $collection_adaptor = $self->reseq_db->get_CollectionAdaptor();
  my $sample_adaptor     = $self->reseq_db->get_SampleAdaptor();

  my $sample_id = $run->sample_id;
  my $sample    = $sample_adaptor->fetch_by_dbID($sample_id);

  my $sample_name        = $sample->sample_alias;
  my $sample_name_reg_ex = qr/\Q$sample_name\E/;

  for my $collection_type ( $self->collection_types_to_check() ) {
    my $collection =
      $collection_adaptor->fetch_by_name_and_type( $run->run_source_id,
      $collection_type );

    next unless $collection;
    throw(
"Found non-file collection for $run->run_source_id $collection_type, cannot check paths"
    ) if ( $collection->table_name ne 'file' );

    for my $file ( @{$collection->others} ) {
      next if ( $file->name =~ $sample_name_reg_ex );

      my $name     = $file->name;
      my $new_name = $file->name;
      my $old_sample;

      if ( $name =~ /(NA\d+)/ ) {
        $old_sample = $1;
      }
      elsif ( $name =~ /unidentified/ ) {
        $old_sample = 'unidentified';
      }
      else {
        throw( "Failed to get a sample id from " . $file->name );
      }

      $new_name =~ s/$old_sample/$sample_name/;
      $new_name =~ s/withdrawn/ftp/;

      if ( $new_name =~ /withdrawn/ ) {
        throw( "Still have withdrawn in new name " . $new_name );
      }

      if ( $name eq $new_name ) {
        throw("There is a problem "
            . $name
            . " is the same as "
            . $new_name
            . " when "
            . "changing from "
            . $old_sample . " to "
            . $sample_name );
      }

      $self->file_hash( $name, $new_name );

    }

  }
}

sub report {
  my ($self) = @_;

  my $file_hash    = $self->file_hash();
  my $sample_count = keys(%$file_hash);
  my $log_fh       = $self->log_fh();

# Delierate use of STDOUT instead of $log_fh - the log_fh should contain the file list, warning about the number of errors should go to stdout
  print STDOUT "There are $sample_count sample issues to resolve$/"
    if ($sample_count);

  while ( my ( $original_path, $suggested_path ) = each %$file_hash ) {
    print $log_fh "$original_path\t$suggested_path$/";
  }
}

sub file_hash {
  my ( $self, $k, $v ) = @_;

  if ( $k && $v ) {
    $self->{file_hash}->{$k} = $v;
  }

  return $self->{file_hash};
}

sub collection_types_to_check {
  return qw(FILTERED_FASTQ);
}

1;
