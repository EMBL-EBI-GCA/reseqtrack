package ReseqTrack::DBSQL::AlignmentMetaInfoAdaptor;

use strict;
use warnings;
use vars qw(@ISA);

use ReseqTrack::DBSQL::BaseAdaptor;
use ReseqTrack::AlignmentMetaInfo;
use ReseqTrack::Tools::Exception qw(throw warning);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::AlignmentMetaInfoUtils qw (are_alignment_meta_infos_identical);

use File::Basename;

@ISA = qw(ReseqTrack::DBSQL::BaseAdaptor);

sub new {
  my ($class, $db) = @_;

  my $self = $class->SUPER::new($db);

  return $self;
}

sub columns{
  return "alignment_meta_info.alignment_meta_info_id, alignment_meta_info.file_id, ".
      "alignment_meta_info.index_file_id, alignment_meta_info.sample_name, ".
      "alignment_meta_info.region, alignment_meta_info.assembly, ".
      "alignment_meta_info.program, alignment_meta_info.mapped_basecount, ".
      "alignment_meta_info.technology";
}

sub table_name{
  return "alignment_meta_info";
}

sub fetch_by_file_name{
  my ($self, $filename) = @_;
  my $sql = "select ".$self->columns." from ".$self->table_name.", file ".
      "where alignment_meta_info.file_id = file.file_id and file.name = ?";
  my $sth = $self->prepare($sql);
  $sth->bind_param(1, $filename);
  $sth->execute;
  my $rowhashref = $sth->fetchrow_hashref;
  my $object = $self->object_from_hashref($rowhashref) if($rowhashref);
  $sth->finish;
  return $object;
}

sub store{
  my ($self, $alignment_meta_info, $update) = @_;
  my $existing = $self->fetch_by_file_name($alignment_meta_info->file->name);
  if($existing){
    $alignment_meta_info->dbID($existing->dbID);
    $alignment_meta_info->adaptor($self);
    if(are_alignment_meta_infos_identical($existing, $alignment_meta_info)){
      $self->store_history($alignment_meta_info);
      $self->store_statistics($alignment_meta_info, $update);
      return $alignment_meta_info;
    }
    if($update){
      return $self->update($alignment_meta_info);
    }
    warning("Failed to store ".$alignment_meta_info." an line pointing to ".
            $alignment_meta_info->file->name." already exists in the database");
    return undef;
  }
  my $sql = "insert into alignment_meta_info(file_id, index_file_id, sample_name, ".
      "region, assembly, program, mapped_basecount, technology) ".
      "values(?, ?, ?, ?, ?, ?, ?, ?)";
  my $sth = $self->prepare($sql);
  #Need to first check if file and index file are already stored
  my $fa = $self->db->get_FileAdaptor;
  unless($alignment_meta_info->file->dbID){
    my $stored = $fa->store($alignment_meta_info->file);
    $alignment_meta_info->file($stored);
  }
  unless($alignment_meta_info->index_file->dbID){
    my $stored = $fa->store($alignment_meta_info->index_file);
    $alignment_meta_info->index_file($stored);
  }
  $sth->bind_param(1, $alignment_meta_info->file_id);
  $sth->bind_param(2, $alignment_meta_info->index_file_id);
  $sth->bind_param(3, $alignment_meta_info->sample_name);
  $sth->bind_param(4, $alignment_meta_info->region);
  $sth->bind_param(5, $alignment_meta_info->assembly);
  $sth->bind_param(6, $alignment_meta_info->program);
  $sth->bind_param(7, $alignment_meta_info->mapped_basecount);
  $sth->bind_param(8, $alignment_meta_info->technology);
  $sth->execute;
  my $dbID = $sth->{'mysql_insertid'};
  $alignment_meta_info->dbID($dbID);
  $alignment_meta_info->adaptor($self);
  $sth->finish;
  $self->store_history($alignment_meta_info);
  $self->store_statistics($alignment_meta_info);
  return $alignment_meta_info;
}

sub object_from_hashref{
  my ($self, $hashref) = @_;
  throw("Can't create object from empty hashref") unless($hashref && keys(%$hashref) >= 1);
  my $alignment_meta_info = ReseqTrack::AlignmentMetaInfo->new(
    -dbID => $hashref->{alignment_meta_info_id},
    -adaptor => $self,
    -file_id => $hashref->{file_id},
    -index_file_id => $hashref->{index_file_id},
    -sample_name => $hashref->{sample_name},
    -region => $hashref->{region},
    -assembly => $hashref->{assembly},
    -program => $hashref->{program},
    -mapped_basecount => $hashref->{mapped_basecount}
      );
  return $alignment_meta_info;
}

1;
