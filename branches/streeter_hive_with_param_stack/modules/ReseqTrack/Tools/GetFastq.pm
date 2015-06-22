=pod

=head1 NAME

ReseqTrack::Tools::GetFastq

=head1 SYNOPSIS

This is a object to get fastq files associated with a run_id

Subroutines in this package are for publicly available fastq files in the ERA database.
Child classes can override any of these subroutines for other methods of getting fastq files.

=head1 Example


=cut

package ReseqTrack::Tools::GetFastq;

use strict;
use warnings;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::ERAUtils qw ();
use ReseqTrack::Tools::FileSystemUtils qw (check_directory_exists delete_file run_md5);
use File::Copy qw (copy);

=head2 new

  Arg [-output_dir]   :
      string, destination directory for fastq files
  Arg [-run_meta_info]   :
      A RunMetaInfo object
  Arg [-db]   :
      A ERADBAdaptor object
  Arg [-source_root_dir]   :
      string, the root directory for era files, e.g. /nfs/era-pub
  Arg [-clobber]   :
      boolean, specifies whether existing files should be copied over or not

=cut


sub new {
  my ($class, @args) = @_;
  my $self ={};
  bless $self,$class;
  my ($output_dir, $run_meta_info, $db, $source_root_dir, $clobber)
    = rearrange([qw(OUTPUT_DIR RUN_META_INFO DB SOURCE_ROOT_DIR CLOBBER
		    )],  @args);

  $self->output_dir($output_dir);
  $self->run_meta_info($run_meta_info);
  $self->db($db);
  $self->source_root_dir($source_root_dir || '/nfs/era-pub');
  $self->clobber($clobber);

  return $self;
}

=head2 run

returns the number of files that were copied

=cut

sub run {
  my $self = shift;

  throw("do not have output_dir") if !$self->output_dir;
  throw("do not have RunMetaInfo object") if !$self->run_meta_info;
  throw("do not have a DB adaptor") if !$self->db;
  throw("do not have a source root directory") if !$self->source_root_dir;

  if ($self->check_fastq_available == 0) {
    throw("no fastq files available for ".$self->run_meta_info->run_id."\n");
  }

  if ($self->check_status == 0) {
    return 0;
  }
  $self->get_fastq_details;
  $self->make_output_hash;
  $self->prepare_output_dir;
  $self->db->dbc->disconnect_when_inactive(1);
  $self->get_files;
  $self->check_sizes;
  $self->check_md5s;
  return scalar @{$self->output_files};
}

=head2 run

returns 1 if fastq are available; otherwise 0

=cut

sub check_fastq_available {
  my $self = shift;
  my $era_rmia = $self->db->get_ERARunMetaInfoAdaptor;
  return $era_rmia->is_fastq_available($self->run_meta_info->run_id);
}

=head2 run

returns 1 if status is OK; otherwise 0

=cut

sub check_status {
  my $self = shift;
  my $era_rmia = $self->db->get_ERARunMetaInfoAdaptor;
  my $status = $era_rmia->get_status($self->run_meta_info->run_id);
  if ($status ne 'public') {
    print "status is not public for ".$self->run_meta_info->run_id."\n";
    return 0;
  }
  return 1;
}

sub get_fastq_details {
  my $self = shift;
  my ($db_md5_hash, $db_size_hash, $name_hash) =
      ReseqTrack::Tools::ERAUtils::get_fastq_details($self->run_meta_info->run_id, $self->db, $self->source_root_dir);
  $self->db_md5_hash($db_md5_hash);
  $self->db_size_hash($db_size_hash);
  $self->name_hash($name_hash);
}

sub make_output_hash {
  my $self = shift;
  my $output_dir = $self->output_dir;
  my $name_hash = $self->name_hash;
  my %output_hash;
  while (my ($basename, $source_path) = each %$name_hash) {
    my $output_path = $output_dir . '/' . $basename;
    $output_path =~ s{//}{/}g;
    $output_hash{$source_path} = $output_path;
  }
  $self->output_hash(\%output_hash);
}

sub prepare_output_dir {
  my $self = shift;
  my $output_dir = $self->output_dir;
  check_directory_exists($output_dir);
  my @exists = grep {-e $_} values %{$self->output_hash};
  if (@exists) {
    if ($self->clobber) {
      print "Deleting files associated with ".$self->run_meta_info->run_id." in $output_dir\n";
      delete_file($_, 1) foreach @exists;
    } else {
      throw("Files already exist and clobber flag is 0: @exists")
    }
  }
}

sub get_files {
  my $self = shift;
  my $output_hash = $self->output_hash;
  while (my($source_path, $output_path) = each %$output_hash) {
    copy($source_path, $output_path)
        or throw("Failed copy from $source_path to $output_path");
    chmod 0644, $output_path;
    $self->output_files($output_path);
  }
}

sub check_md5s {
  my $self = shift;
  my $md5_hash = $self->db_md5_hash;
  my $output_hash = $self->output_hash;
  while (my($source_path, $output_path) = each %$output_hash) {
    my $output_md5 = run_md5($output_path);
    throw("Problem with output md5 for $source_path: the database gives $md5_hash->{$source_path} but run_md5 gives $output_md5")
        if ($output_md5 ne $md5_hash->{$source_path});
    $self->md5_hash->{$output_path} = $output_md5;
  }
}

sub check_sizes {
  my $self = shift;
  my $size_hash = $self->db_size_hash;
  my $output_hash = $self->output_hash;
  while (my($source_path, $output_path) = each %$output_hash) {
    my $output_size = -s $output_path;
    throw("Problem with output size for $source_path: the database gives $size_hash->{$source_path} but perl gives $output_size")
        if ($output_size != $size_hash->{$source_path});
    $self->size_hash->{$output_path} = $output_size;
  }
}





=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : Various things, ints, strings, objects
  Function  : hold object and return on request
  Returntype: Various things
  Exceptions: some methods throw if not given the correct object/type or
  if the variable doesn't match some specified criteria like a non existent file
  or directory
  Example   : my $output_dir = $self->output_dir

=cut

sub output_dir{
  my ($self, $output_dir) = @_;
  if($output_dir){
    $self->{output_dir} = $output_dir;
  }
  $self->{output_dir};
}

sub run_meta_info{
  my ($self, $run_meta_info) = @_;
  if($run_meta_info){
    $self->{run_meta_info} = $run_meta_info;
  }
  return $self->{run_meta_info};
}

sub db{
  my ($self, $db) = @_;
  if($db){
    $self->{db} = $db;
  }
  return $self->{db};
}

sub source_root_dir{
  my ($self, $source_root_dir) = @_;
  if($source_root_dir){
    $self->{source_root_dir} = $source_root_dir;
  }
  return $self->{source_root_dir};
}

sub db_md5_hash{
  my ($self, $db_md5_hash) = @_;
  $self->{db_md5_hash} ||= {};
  if($db_md5_hash){
    $self->{db_md5_hash} = $db_md5_hash;
  }
  return $self->{db_md5_hash};
}

sub md5_hash{
  my ($self, $md5_hash) = @_;
  $self->{md5_hash} ||= {};
  if($md5_hash){
    $self->{md5_hash} = $md5_hash;
  }
  return $self->{md5_hash};
}

sub size_hash{
  my ($self, $size_hash) = @_;
  $self->{size_hash} ||= {};
  if($size_hash){
    $self->{size_hash} = $size_hash;
  }
  return $self->{size_hash};
}

sub db_size_hash{
  my ($self, $db_size_hash) = @_;
  $self->{db_size_hash} ||= {};
  if($db_size_hash){
    $self->{db_size_hash} = $db_size_hash;
  }
  return $self->{db_size_hash};
}


sub name_hash{
  my ($self, $name_hash) = @_;
  $self->{name_hash} ||= {};
  if($name_hash){
    $self->{name_hash} = $name_hash;
  }
  return $self->{name_hash};
}

sub output_hash{
  my ($self, $output_hash) = @_;
  $self->{output_hash} ||= {};
  if($output_hash){
    $self->{output_hash} = $output_hash;
  }
  return $self->{output_hash};
}

sub clobber{
  my ($self, $clobber) = @_;
  if(defined($clobber)){
    $self->{clobber} = $clobber;
  }
  return $self->{clobber};
}

sub output_files {
  my $self = shift;
  $self->{output_files} ||= [];
  if (@_) {
    push(@{$self->{output_files}}, grep {$_} @_);
  }
  return $self->{output_files};
}

1;
