package ReseqTrack::Tools::AlignmentMetaInfoUtils;

use strict;
use warnings;
use Exporter;
use ReseqTrack::AlignmentMetaInfo;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileUtils;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(are_alignment_meta_infos_identical);



=head2 are_alignment_meta_infos_identical

  Arg [1]   : ReseqTrack::AlignmentMetaInfo
  Arg [2]   : ReseqTrack::AlignmentMetaInfo
  Function  : compare two objects return 0 if they aren't the same and 1 if they are
  the only aspects which aren't compared are updated and history
  Returntype: 0/1
  Exceptions: none
  Example   : if(are_alignment_meta_infos_identical($file, $other_file){
     print $file." is the same as ".$other_file."\n";
   }

=cut


sub are_alignment_meta_infos_identical{
  my ($one, $two) = @_;

  throw("Must pass are_alignment_meta_infos_identical two AlignmentMetaInfo objects") unless($one && $two);
  throw("Must pass are_alignment_meta_infos_identical two AlignmentMetaInfo ".
        "objects and not ".$one." and ".$two) 
      unless($one->isa("ReseqTrack::AlignmentMetaInfo") && 
             $two->isa("ReseqTrack::AlignmentMetaInfo"));
  if($one->file && $two->file){
    return 0 unless(are_files_identical($one->file, $two->file));
  }
  return 0 unless($one->file_id == $two->file_id);
  return 0 unless($one->index_file_id == $two->index_file_id);
  if($one->index_file && $two->index_file){
    return 0 unless(are_files_identical($one->index_file, $two->index_file));
  }
  return 0 unless($one->sample_name eq $two->sample_name);
  return 0 unless($one->region eq $two->region);
  return 0 unless($one->assembly eq $two->assembly);
  return 0 unless($one->program eq $two->program);
  return 0 unless($one->technology eq $two->technology);
  return 0 unless($one->mapped_basecount == $two->mapped_basecount);
                  
  return 1;
}








1;
