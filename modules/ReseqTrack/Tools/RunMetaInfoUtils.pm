package ReseqTrack::Tools::RunMetaInfoUtils;

use strict;
use Exporter;
use ReseqTrack::RunMetaInfo;
use ReseqTrack::History;
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::GeneralUtils qw(current_time);
use File::Basename;

use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(are_run_meta_infos_identical create_object_from_index_line
        index_to_name_hash get_meta_info_hash create_history_for_run_meta_info
        create_suppressed_index_line create_index_line
        get_files_associated_with_run get_file_collections_associated_with_run
        copy_run_meta_info get_analysis_group get_sequence_index_stats
        get_withdrawn_summary get_study_descriptions get_index_group_stats
        get_run_id_from_filename convert_center_name convert_population
        create_directory_path);



=head2 are_run_meta_infos_identical

  Arg [1]   : ReseqTrack::RunMetaInfo
  Arg [2]   : ReseqTrack::RunMetaInfo
  Function  : compare two objects return 0 if they aren't the same and 1 if they are
  the only aspects which aren't compared are updated and history
  Returntype: 0/1
  Exceptions: none
  Example   : if(are_run_meta_infos_identical($file, $other_file){
     print $file." is the same as ".$other_file."\n";
   }

=cut


sub are_run_meta_infos_identical{
  my ($one, $two, $skip_date) = @_;
  my $verbose = 0;
  throw("Must pass are_run_meta_infos_identical two RunMetaInfo objects") unless($one && $two);
  throw("Must pass are_run_meta_infos_identical two RunMetaInfo objects and not ".$one." and ".
        $two) unless($one->isa("ReseqTrack::RunMetaInfo") && 
                     $two->isa("ReseqTrack::RunMetaInfo"));
  
  
  return 0 if($one->run_id ne $two->run_id);
  return 0 if($one->study_id ne $two->study_id);
  return 0 if($one->study_name ne $two->study_name);
  return 0 if(($one->center_name || $two->center_name) && 
              ($one->center_name ne $two->center_name));
  return 0 if($one->submission_id ne $two->submission_id);
  return 0 if($one->sample_id ne $two->sample_id);
  return 0 if($one->library_strategy ne $two->library_strategy);
  unless($skip_date){
    if($two->submission_date){
      if($one->submission_date){
        if($one->submission_date ne $two->submission_date){
          return 0;
        }
      }
    }
  }
  if(!$one->sample_name || !$two->sample_name){
    print STDERR $one->run_id." or ".$two->run_id." doesn't have a sample defined\n";
  }
  return 0 if($one->sample_name ne $two->sample_name);
  if(($one->population || $two->population) && 
     ($one->population ne $two->population)){
    if($one->population){
      return 0;
    }
  }
  return 0 if($one->instrument_platform ne $two->instrument_platform);
  return 0 if($one->experiment_id ne $two->experiment_id);
  if(($one->instrument_model && $two->instrument_model) && ($one->instrument_model ne $two->instrument_model)){
    return 0;
  }
  return 0 if($one->library_name ne $two->library_name);
  return 0 if($one->run_name ne $two->run_name);
  if($one->run_block_name || $two->run_block_name){
    unless($one->run_block_name eq "NULL"){
      return 0 unless($one->run_block_name eq $two->run_block_name);
    }
  }
  if(($one->status || $two->status) &&
     ($one->status ne $two->status)){
    return 0;
  }
  return 0 if($one->paired_length != $two->paired_length);
  return 0 if($one->library_layout ne $two->library_layout);

  if($one->archive_base_count || $two->archive_base_count){
    if($one->archive_base_count ne $two->archive_base_count){
      unless(!$one->archive_base_count || $one->archive_base_count == 0){
        print "ONE ".$one->archive_base_count." is different to TWO ".$two->archive_base_count.
            " returning 0\n" if($verbose);
        return 0;
      }
    }
  }
  if($one->archive_read_count || $two->archive_read_count){
    if($one->archive_read_count ne $two->archive_read_count){
      return 0 unless(!$one->archive_read_count || $one->archive_read_count == 0);
    }
  }
 
  return 1;

}

sub create_history_for_run_meta_info{
  my ($new, $old) = @_;
  throw("Must pass create_history_for_run_meta_info two RunMetaInfo objects") unless($new && $old);
  throw("Must pass create_history_for_run_meta_info two RunMetaInfo objects and not ".$new." and ".
        $old) unless($new->isa("ReseqTrack::RunMetaInfo") && 
                     $old->isa("ReseqTrack::RunMetaInfo"));
  
  
  throw("Can't create a history object if you are changing the run id from ".
        $old->run_id." to ".$new->run_id) unless($new->run_id eq $old->run_id);
  my $comment;
  if($old->submission_date && $new->submission_date){
    $comment = "Changing submission date from ".$old->submission_date." to ".$new->submission_date if($old->submission_date ne $new->submission_date);
  }
  $comment = "Changing library strategy ".$old->library_strategy." to ".$new->library_strategy if($old->library_strategy ne $new->library_strategy);
    $comment = "Changing study name ".$old->study_name." ".$new->study_name if($new->study_name ne $old->study_name);
  
  $comment = "Changing submission id ".$old->submission_id." ".$new->submission_id if($new->submission_id ne $old->submission_id);
  
  
  $comment = "Changing population ".$old->population." ".$new->population 
      if(($new->population || $old->population) && 
         ($new->population ne $old->population));
  $comment = "Changing instrument platform ".$old->instrument_platform." ".
      $new->instrument_platform 
      if($new->instrument_platform ne $old->instrument_platform);
  $comment = "Changing experiment id ".$old->experiment_id." ".$new->experiment_id
      if($new->experiment_id ne $old->experiment_id);
  if(($new->instrument_model && $old->instrument_model) 
     && ($new->instrument_model ne $old->instrument_model)){
    $comment = "Changing instrument model ".$old->instrument_model." ".
        $new->instrument_model;
  }
  $comment = "Changing library name ".$old->library_name." ".$new->library_name if($new->library_name ne $old->library_name);
  $comment = "Changing run name ".$old->run_name." ".$new->run_name if($new->run_name ne $old->run_name);
  if($old->run_block_name){
    if(($new->run_block_name || $old->run_block_name) && 
       ($new->run_block_name ne $old->run_block_name)){
      $comment = "Changing run block name ".$old->run_block_name." ".$new->run_block_name;
    }
  }
  if(($new->status || $old->status) &&
     ($new->status ne $old->status)){
    $comment = "Changing status ".$old->status." ".$new->status;
  }
  $comment = "Changing paired length ".$new->paired_length." ".$old->paired_length  if($new->paired_length != $old->paired_length);
  $comment = "Changing library layout ".$old->library_layout." ".$new->library_layout if($new->library_layout ne $old->library_layout);
  $comment = "Changing archive base count from ".$old->archive_base_count." to ".$new->archive_base_count if($new->archive_base_count != $old->archive_base_count);
  $comment = "Changing archive read count from ".$old->archive_read_count." to ".$new->archive_read_count if($new->archive_read_count != $old->archive_read_count);
  $comment = "Changing study id ".$old->study_id." ".$new->study_id if($new->study_id ne $old->study_id);
  $comment = "Changing center_name ".$old->center_name." ".$new->center_name 
      if(($new->center_name || $old->center_name) && 
         ($new->center_name ne $old->center_name));
  $comment = "Changing sample id ".$old->sample_id." ".$new->sample_id  
      if($new->sample_id ne $old->sample_id);
  if(!$new->sample_name || !$old->sample_name){
    print STDERR $new->run_id." or ".$old->run_id." doesn't have a sample defined\n";
  }
  $comment = "Changing sample name ".$old->sample_name." ".$new->sample_name if($new->sample_name ne $old->sample_name);
  if(!$comment){
    #print "There are no differences between ".$old." and ".$new." can't create a history\n";
    return undef;
  }
  my $history = ReseqTrack::History->new(
    -other_id => $new->dbID,
    -table_name => 'run_meta_info',
    -comment => $comment,
      );
  return $history;
}


sub create_object_from_index_line{
  my ($line) = @_;
  my @values = split /\t/, $line;
  my $run_meta_info =  ReseqTrack::RunMetaInfo->new(
    -run_id => $values[2],
    -study_id => $values[3],
    -study_name => $values[4],
    -center_name => $values[5],
    -submission_id => $values[6],
    -submission_date => $values[7],
    -sample_id => $values[8],
    -sample_name => $values[9],
    -population => $values[10],
    -experiment_id => $values[11],
    -instrument_platform => $values[12],
    -instrument_model => $values[13],
    -library_name => $values[14],
    -run_name => $values[15],
    -run_block_name => $values[16],
    -library_layout => $values[18],
    -paired_length => $values[17],
      );
  return $run_meta_info;
}

sub copy_run_meta_info{
  my ($rmi) = @_;
  my $run_meta_info =  ReseqTrack::RunMetaInfo->new(
    -run_id => $rmi->run_id,
    -study_id => $rmi->study_id,
    -study_name => $rmi->study_name,
    -center_name => $rmi->center_name,
    -submission_id => $rmi->submission_id,
    -submission_date => $rmi->submission_date,
    -sample_id => $rmi->sample_id,
    -sample_name => $rmi->sample_name,
    -population => $rmi->population,
    -experiment_id => $rmi->experiment_id,
    -instrument_platform => $rmi->instrument_platform,
    -instrument_model => $rmi->instrument_model,
    -library_name => $rmi->library_name,
    -run_name => $rmi->run_name,
    -run_block_name => $rmi->run_block_name,
    -library_layout => $rmi->library_layout,
    -paired_length => $rmi->paired_length, 
    -status => $rmi->status,
      );
  return $run_meta_info;
}

sub index_to_name_hash{
  my %hash;
  $hash{0} = 'fastq_file';
  $hash{1} = 'md5';
  $hash{2} = 'run_id';
  $hash{3} = 'study_id';
  $hash{4} = 'study_name';
  $hash{5} = 'center_name';
  $hash{6} = 'submission_id';
  $hash{7} = 'submission_date';
  $hash{8} = 'sample_id';
  $hash{9} = 'sample_name';
  $hash{10} = 'population';
  $hash{11} = 'experiment_id';
  $hash{12} = 'instrument_platform';
  $hash{13} = 'instrument_model';
  $hash{14} = 'library_name';
  $hash{15} = 'run_name';
  $hash{16} = 'run_block_name';
  $hash{17} = 'paired_length';
  $hash{18} = 'library_layout';
  $hash{19} = 'paired_fastq';
  $hash{20} = 'withdrawn';
  $hash{21} = 'withdrawn_date';
  $hash{22} = 'comment';
  $hash{23} = 'archive_read_count';
  $hash{24} = 'archive_base_count';
  return \%hash;
}


sub get_meta_info_hash{
  my ($meta_info) = @_;
  my %hash;
  foreach my $object(@$meta_info){
    next unless($object);
    if($hash{$object->run_id}){
      #warning($object->run_id." already exists in hash skipping");
      next;
    }
    $hash{$object->run_id} = $object;
  }
  return \%hash;
}

sub index_method_array{
  return ['run_id', 'study_id', 'study_name', 'center_name', 'submission_id', 
          'submission_date', 'sample_id', 'sample_name', 'population', 
          'experiment_id', 'instrument_platform', 'instrument_model', 
          'library_name', 'run_name', 'run_block_name', 'paired_length',
          'library_layout'];
}

sub create_index_line{
  my ($file, $md5, $rmi, $paired_file, $withdrawn, $comment, $time, $read_count, $base_count, $analysis_group) = @_;
  no warnings;
  if($withdrawn){
    $time = current_time() if(!$time);
  }
  unless($analysis_group){
    $analysis_group = get_analysis_group($rmi);
  }
  $withdrawn = 0 unless($withdrawn);
  $read_count = $rmi->archive_read_count if(!$read_count);
  $read_count = 'not available' if (!$read_count);
  $base_count = $rmi->archive_base_count if(!$base_count);
  $base_count = 'not available' if (!$base_count);
  my $methods = index_method_array();
  my $string;
  $string .=  $file."\t";
  $string .= $md5."\t";
  foreach my $method(@$methods){
    $string .= $rmi->$method."\t";
  }
  $string .= $paired_file."\t";
  $string .= $withdrawn."\t";
  $string .= $time."\t";
  $string .= $comment."\t";
  $string .= $read_count."\t";
  #$string .= $base_count."\t";
  $string .= $base_count."\t";
  $string .= $analysis_group;
  return $string;
}

sub create_suppressed_index_line{
  my ($rmi, $comment,$time, $analysis_group) = @_;
  unless($rmi->status ne 'public'){
    #warning($rmi->run_id." has status ".$rmi->status." are you sure you want to ".
    #    "print a suppressed style line for it");
  }
  $comment = 'SUPPRESSED IN ARCHIVE' unless($comment);
  my $md5 = "................................";
  my $file = "data/".$rmi->sample_name."/sequence_read/".$rmi->run_id.".fastq.gz";

  return create_index_line($file, $md5, $rmi, undef, 1, $comment, $time, undef, 
                           undef, $analysis_group);
}

sub get_files_associated_with_run{
  my ($rmi, $type) = @_;
  throw("Can't fetch files associated with ".$rmi->run_id." it has no adaptor")
      unless($rmi->adaptor);
  my $fa = $rmi->adaptor->db->get_FileAdaptor;
  my $files = $fa->fetch_all_like_name($rmi->run_id);
  my @return_files;
  if($type){
    foreach my $file(@$files){
      next if($file->name =~ /\/withdrawn\//);
      if($file->type eq $type){
        push(@return_files, $file);
      }
    }
    return \@return_files;
  }else{
    return $files;
  }
}

sub get_file_collections_associated_with_run{
  my ($rmi, $type) = @_;
  throw("Can't fetch collections associated with ".$rmi->run_id." it has no adaptor")
      unless($rmi->adaptor);
  my $fa = $rmi->adaptor->db->get_CollectionAdaptor;
  my $collections = $fa->fetch_by_name_and_table_name($rmi->run_id, 'file');
  my @return_collections;
  if($type){
    foreach my $collection(@$collections){
      if($collection->type eq $type){
        push(@return_collections, $collection);
      }
    }
    return \@return_collections;
  }else{
    return $collections;
  }
}


sub get_analysis_group{
  my ($rmi) = @_;
  my $analysis_group = 'low coverage';
  if($rmi->study_id eq 'SRP000032'){
    $analysis_group = 'high coverage';
    }
  if($rmi->study_id eq 'SRP000033'){
    $analysis_group = 'exon targetted';
  }
  return $analysis_group;
}

sub get_sequence_index_stats{
  my ($index_file, $summary_column, $stats_column, $run_id_hash, $population_rules) = @_;
  open(FH, $index_file) or 
      throw("ReseqTrack::Tools::RunMetaInfoUtils::get_sequence_index_stats ".
            "failed to open ".$index_file." $!");
  my %hash;
  while(<FH>){
    next if(/SUBMISSION_ID/);
    chomp;
    my @values = split /\t/;
    next if($run_id_hash && !($run_id_hash->{$values[2]}));
    $values[10] = convert_population($values[10], $population_rules, $values[2], undef);
    $values[5] = convert_center_name($values[5]);
    next if($values[20]);
    my $summary = $values[$summary_column];
    #print "study is $summary\trun_id is $values[2]\n";
    my $stat = $values[$stats_column];
    throw("Can't get summary stats from ".$index_file." based on ".$summary_column.
          " as ".$stats_column." contains some non numeric values") 
        unless($stat =~ /^\d+$/);
    $hash{$summary} += $stat;
  }
  close(FH);
  return \%hash;
}

sub get_index_group_stats{
  my ($index_file, $first_column, $second_column, $stats_column, $run_id_hash, $population_rules) = @_;
   open(FH, $index_file) or 
      throw("ReseqTrack::Tools::RunMetaInfoUtils::get_sequence_index_stats ".
            "failed to open ".$index_file." $!");
  my %hash;
  while(<FH>){
    next if(/SUBMISSION_ID/);
    chomp;
    my @values = split /\t/;
    next if($run_id_hash && !($run_id_hash->{$values[2]}));
    $values[10] = convert_population($values[10], $population_rules, $values[2], undef);
    $values[5] = convert_center_name($values[5]);
    next if($values[20]);
    $hash{$values[$first_column]}{$values[$second_column]} += $values[$stats_column];
  }
  close(FH);
  return \%hash;
}

sub get_withdrawn_summary{
  my ($index_file, $run_id_hash) = @_;
  open(FH, $index_file) or 
      throw("ReseqTrack::Tools::RunMetaInfoUtils::get_withdrawn_summary ".
            "failed to open ".$index_file." $!");
  my %hash;
  while(<FH>){
    next if(/SUBMISSION_ID/);
    chomp;
    my @values = split /\t/;
    next if($run_id_hash && !($run_id_hash->{$values[2]}));
    next unless($values[20]);
    my $reason = $values[22];
    $hash{$reason}++;
  }
  close(FH);
  return \%hash;
}

sub get_study_descriptions{
  my ($index_file) = @_;
  open(FH, $index_file) or 
      throw("ReseqTrack::Tools::RunMetaInfoUtils::get_study_descriptions ".
            "failed to open ".$index_file." $!");
  my %hash;
  while(<FH>){
    next if(/SUBMISSION_ID/);
    chomp;
    my @values = split /\t/;
    $hash{$values[3]} = $values[4];
  }
  close(FH);
  return \%hash;
}


sub convert_center_name{
  my ($name) = @_;
  if($name && $name eq 'ABHTD'){
    return 'ABI';
  }else{
    return $name;
  }
}



sub convert_population{
  my ($string, $population_rules, $run_id, $study_id) = @_;
  throw("Can't convert an empty string for ".$run_id." ".$study_id)
    unless($string);
  foreach my $pr (@$population_rules) {
    if ($pr->test_description($string)) {
      return $pr->population;
    }
  }
  throw("Failed to find pop for ".$string." ".$run_id." ".$study_id);
}


sub get_run_id_from_filename{
  my ($name) = @_;
  my $filename = basename($name);
  $filename =~ /([E|S]RR\d+)/;
  return $1;
}

=head2 create_directory_path

  Arg [1]   : ReseqTrack::RunMetaInfo
  Arg [2]   : string, directory layout e.g. 'sample_name/archive_sequence', tokens matching method names in RunMetaInfo will be substituted with that method's return value
  Arg [3]   : string, base directory
  Function  : Creates a directory path, combining the base directory path and run_meta_info values
  Returntype: string
  Exceptions: none
  Example   : $dir = create_directory_path($rmi, 'population/sample_name', '/path/to/dir')

=cut

sub create_directory_path{
  my ($run_meta_info, $directory_layout, $base_directory) = @_;
  my $dir_path = $base_directory;
  my @layout_chunks =  split( /\//, $directory_layout);
  my $method_matches = 0;
  foreach my $layout_chunk (@layout_chunks) {
    my $method = $run_meta_info->can($layout_chunk);	
    if ($method){
      $layout_chunk = &$method($run_meta_info);
      $method_matches++;
    }
            
    $dir_path .= '/' if ($dir_path);
    $dir_path .= $layout_chunk;
  }

  throw "Directory layout ($directory_layout) did not call any run_meta_info methods"
    if (@layout_chunks && ! $method_matches);

  $dir_path =~ s/\/\//\//;
  return $dir_path;
}


1;
