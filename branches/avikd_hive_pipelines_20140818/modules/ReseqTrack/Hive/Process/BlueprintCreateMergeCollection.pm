package ReseqTrack::Hive::Process::BlueprintCreateMergeCollection;

use strict;

use base ('ReseqTrack::Hive::Process::BaseProcess');
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception qw(throw);

use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_directory_exists check_file_does_not_exist);
use ReseqTrack::Collection;
use File::Basename qw( dirname );
use List::Compare;
use IPC::System::Simple qw(capture);

sub run {
    my $self = shift @_;

    $self->param_required('input_bam');
    my $input_bams = $self->param_as_array('input_bam');
    throw('no input bams') 
         if !@$input_bams;

    $self->param_required('bam');  
    my $merged_bams = $self->param_as_array('bam');
    throw('got multiple merged bams per experiment')   
         if scalar @$merged_bams > 1;

    my $samtools = $self->param_required('samtools');
    my $collection_type = $self->param_required('collection_type');
    my $collection_name = $self->param_required('collection_name');
    my $tag_name = $self->param_required('rg_tag_name');
    my $db_params = $self->param_required('reseqtrack_db');
    my $db = ReseqTrack::DBSQL::DBAdaptor->new(%{$db_params});
    my $fa = $db->get_FileAdaptor;
    my $ca = $db->get_CollectionAdaptor;
  
    my $merged_bam = $$merged_bams[0];
    check_file_exists( $merged_bam );
    my $merged_bam_rg_ids = _get_rg_tags( $merged_bam, $tag_name, $samtools ); ## list of runs in the exp bam
    my %run_rg_ids_hash;

    my @new_file_objects;
    my @existing_file_objects;

    foreach my $input_bam ( @$input_bams ){
      check_file_exists( $input_bam );
      my $file = $fa->fetch_by_name( $input_bam );
      my $file_id = $file->dbID;
     
      my $file_collections;
      $file_collections = $ca->fetch_by_other_id_and_type( $file_id, $collection_type)
          if ( $ca->fetch_by_other_id_and_type( $file_id, $collection_type));

      if ( scalar @$file_collections > 0 ){  ## existing merge collection for the same files
        push( @existing_file_objects, $file);             
      }
      else {
        push( @new_file_objects, $file );
      } 

      my $run_bam_rg_ids = _get_rg_tags( $input_bam, $tag_name, $samtools );
      %run_rg_ids_hash = map { $_ => 1 } @$run_bam_rg_ids; 
    }  
    
    my @run_rg_ids_list = map { $_ } sort { $a<=>$b } keys %run_rg_ids_hash;


    my $list_comp = List::Compare->new( {
                         lists    => [ $merged_bam_rg_ids, \@run_rg_ids_list],
                         unsorted => 1,
                    });

    my @merge_bam_only = $list_comp->get_Lonly;
    my @run_bam_only = $list_comp->get_Ronly;
   
    throw('RG ids are not same in the merged bam and run level bams')
             if ( scalar @merge_bam_only > 1 || scalar @run_bam_only > 1 );     ## RG ids shouls be same

    if ( scalar @existing_file_objects == 0 && scalar @new_file_objects >0 ){   ## only update new files
      _create_collection( $db, $ca, \@new_file_objects, $collection_name, $collection_type);  
    }
    elsif ( scalar @existing_file_objects == 0 && scalar @new_file_objects == 0) { ## no files found
      throw('no file found');
    }  
    elsif ( scalar @existing_file_objects > 0) { ## implement methods for existing files collection
      throw("existing collection found, not suported yet: $collection_name & $collection_type :". scalar @existing_file_objects . ":". $$input_bams[0] .":". $$input_bams[1]);
    } 
    $db->dbc->disconnect_when_inactive(1);
}

sub _create_collection {
  my ( $db, $ca, $files, $collection_name, $type ) = @_;
  
  my $collection = ReseqTrack::Collection->new(  -name => $collection_name, 
                                                 -type => $type,
                                                 -others => $files, 
                                                 -table_name =>'file',
                    );  
  $ca->store($collection); 
}

sub _get_rg_tags{
  my ( $bam, $tag_name, $samtools ) = @_;
  my $sam_rg_cmd = "$samtools view -H ";
  $sam_rg_cmd .=  $bam;

  my @header_lines = capture( $sam_rg_cmd );
  my @rg_lines = grep{ /^\@RG/ } @header_lines;

  throw("no RG tag found in $bam")
     unless scalar @rg_lines > 0;

  my @rg_values = split '\s+', $rg_lines[0];
  my @rg_tags = grep{ /^$tag_name/ } @rg_values;
  my @tag_lists;

  foreach my $rg_tag ( @rg_tags ){
    my ($tag, $value) =  split ':', $rg_tag;
    throw("$tag_name not found in the RG tag of $bam")
       unless $value;

    push( @tag_lists, $value);
  }
  return \@tag_lists;
}


1;
