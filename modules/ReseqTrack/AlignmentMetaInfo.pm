=pod

=head1 NAME

ReseqTrack::AlignmentMetaInfo;

=head1 SYNOPSIS

An container object to hold the data for the aligmment_meta_info table for the 
ReseqTrack database

my $alignment_meta_info = ReseqTrack::AlignmentMetaInfo->new(
  -file => $file_object,
  -index_file => $index_file,
  -sample_name => 'NA12878',
  -assembly => 'NCBI36',
  -program => 'maq',
  -technology => 'SLX',
 );

This object is mostly just a container for the database information. 
The functionality it does contain allows it to fetch file objects on the basis of 
given file object id so you can equally create the above object like

my $alignment_meta_info = ReseqTrack::AlignmentMetaInfo->new(
  -file_id => 1,
  -index_file_id => 2,
  -sample_name => 'NA12878',
  -assembly => 'NCBI36',
  -program => 'maq',
  -technology => 'SLX',
  -adaptor => $alignment_meta_info_adaptor
 );

Then when you which to retrive the specified objects the code is capable of using the
specified adaptor to do so

The object inherits from the ReseqTrack::HasHistory object. HasHistory provides some
basic methods for handling both ReseqTrack::History and ReseqTrack::Statistic objects

=cut

package ReseqTrack::AlignmentMetaInfo;

use strict;
use warnings;
use vars qw(@ISA);
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);


use ReseqTrack::HasHistory;


@ISA = qw(ReseqTrack::HasHistory);



=head2 new

  Arg [1]   : ReseqTrack::AlignmentMetaInfo
  Arg [2]   : int, file_id,
  Arg [3]   : ReseqTrack::File
  Arg [4]   : string, sample name, normally in the form NAXXXXXX
  Arg [5]   : string, chromosome number
  Arg [6]   : string, assembly
  Arg [7]   : string, program name
  Arg [8]   : int, number of mapped bases in bam
  Arg [9]   : ReseqTrack::File, for bai file
  Arg [10]  : int, file id for bai file
  Arg [11]  : string, the instrument platform the sequence is from
  Function  : create a ReseqTrack::AlignmentMetaInfo object
  Returntype: ReseqTrack::AlignmentMetaInfo
  Exceptions: 
  Example   : 

=cut



sub new{
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($file_id, $file, $sample_name, $region, $assembly, $program, $mapped_basecount,
      $index_file, $index_file_id, $technology) = 
      rearrange([qw(FILE_ID FILE SAMPLE_NAME REGION ASSEMBLY PROGRAM 
  MAPPED_BASECOUNT INDEX_FILE INDEX_FILE_ID TECHNOLOGY)], @args);
  #Setting defaults
  $region = 'genome' unless($region);
  ######

  $self->file_id($file_id);
  $self->file($file);
  $self->index_file_id($index_file_id);
  $self->index_file($index_file);
  $self->sample_name($sample_name);
  $self->region($region);
  $self->assembly($assembly);
  $self->program($program);
  $self->mapped_basecount($mapped_basecount);
  $self->technology($technology);

  return $self
}



=head2 Accessor methods

  Arg [1]   : ReseqTrack::AlignmentMetaInfo
  Arg [2]   : string/int
  Function  : defined variable 
  Returntype: string/int
  Exceptions: 
  Example   : 

=cut



sub sample_name{
  my ($self, $sample_name) = @_;
  if($sample_name){
    $self->{sample_name} = $sample_name;
  }
  return $self->{sample_name};
}

sub region{
  my ($self, $region) = @_;
  if($region){
    $self->{region} = $region;
  }
  return $self->{region};
}

sub assembly{
  my ($self, $assembly) = @_;
  if($assembly){
    $self->{assembly} = $assembly;
  }
  return $self->{assembly};
}

sub program{
  my ($self, $program) = @_;
  if($program){
    $self->{program} = $program;
  }
  return $self->{program};
}

sub mapped_basecount{
  my ($self, $mapped_basecount) = @_;
  if($mapped_basecount){
    $self->{mapped_basecount} = $mapped_basecount;
  }
  return $self->{mapped_basecount};
}

sub file_id{
  my ($self, $file_id) = @_;
  if($file_id){
    $self->{file_id} = $file_id;
  }
  if(!$self->{file_id} && $self->{file}){
    $self->{file_id} = $self->{file}->dbID;
  }
  return $self->{file_id};

}
sub file{
  my ($self, $file) = @_;
  if($file){
    $self->{file} = $file;
  }
  unless($self->{file}){
    if($self->file_id && $self->adaptor){
      $file = $self->get_file_object($self->file_id);
      $self->{file} = $file;
    }
  }
  return $self->{file};
}
sub index_file_id{
  my ($self, $index_file_id) = @_;
  if($index_file_id){
    $self->{index_file_id} = $index_file_id;
  }
  if(!$self->{index_file_id} && $self->{index_file}){
    $self->{index_file_id} = $self->{index_file}->dbID;
  }
  return $self->{index_file_id};

}
sub index_file{
  my ($self, $index_file) = @_;
  if($index_file){
    $self->{index_file} = $index_file;
  }
  unless($self->{index_file}){
    if($self->index_file_id && $self->adaptor){
      $index_file = $self->get_file_object($self->index_file_id);
      $self->{index_file} = $index_file;
    }
  }
  return $self->{index_file};
}

sub get_file_object{
  my ($self, $dbID) = @_;
  return undef unless($self->adaptor);
  $dbID = $self->file_id unless($dbID);
  my $fa = $self->adaptor->db->get_FileAdaptor;
  return $fa->fetch_by_dbID($dbID);
}

sub technology{
  my ($self, $technology) = @_;
  if($technology){
    $self->{technology} = $technology;
  }
  return $self->{technology};
}


=head2 object_table_name

  Arg [1]   : ReseqTrack::AlignmentMetaInto
  Function  : return table name which represents object, alignment_meta_info
  Returntype: string
  Exceptions: 
  Example   : 

=cut


sub object_table_name{
  my ($self) = @_;
  return 'alignment_meta_info';
}

1;
