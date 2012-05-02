=pod

=head1 NAME

ReseqTrack::Tools::FilterFastq

=head1 SYNOPSIS

This is a object to filter fastq files on the basis of length of read, 
quality values and percent Ns

=head1 Example


=cut

package ReseqTrack::Tools::FilterFastq;

use strict;
use warnings;
use Scalar::Util qw(looks_like_number);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::Scalar qw(assert_ref);
use ReseqTrack::Tools::SequenceIndexUtils qw(assign_files);
use ReseqTrack::Tools::FileSystemUtils;
use FileHandle;
use File::Basename;

=head2 new

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : ReseqTrack::Collection, collection of fastq files
  Arg [3]   : Int, minimum length of sequence
  Arg [4]   : Int, minimum qual value of sequence
  Arg [5]   : Int, maxium percentage of Ns/1 base type in sequence
  Arg [6]   : String, path to output dir
  Arg [7]   : boolean, true if sequence is colorspace
  Arg [8]   : String, run id of sequence must match string in header of each fastq 
              entry
  Arg [9]   : boolean, true if to over write existing files
  Arg [10]  : arrayref of file paths
  Arg [11]  : string, path to mate1
  Arg [12]  : string, path to mate2
  Arg [13]  : string, path to frag file
  Function  : Create a FilterFastq object
  Returntype: ReseqTrack::Tools::FilterFastq
  Exceptions: 
  Example   : 

=cut


sub new {
  my ($class, @args) = @_;
  my $self ={};
  bless $self,$class;
  my ($collection, $min_length, $min_qual, $max_percent_n, $output_dir, 
      $is_colorspace, $run_id, $clobber, $input_files, $mate1, $mate2, $frag) 
    = rearrange([qw(COLLECTION MIN_LENGTH MIN_QUAL MAX_PERCENT_N
		    OUTPUT_DIR IS_COLORSPACE RUN_ID CLOBBER INPUT_FILES
		    MATE1 MATE2 FRAG)],  @args);

  #SETTING DEFAULTS
  $self->min_length(25);
  $self->min_qual(2);
  $self->max_percent_n(0.5);
  $self->clobber(1);
  ########

  $self->collection($collection);
  $self->min_length($min_length);
  $self->min_qual($min_qual);
  $self->max_percent_n($max_percent_n);
  $self->is_colorspace($is_colorspace);
  $self->output_dir($output_dir);
  $self->run_id($run_id);
  $self->clobber($clobber);
  $self->input_files($input_files);
  $self->mate1($mate1);
  $self->mate2($mate2);
  $self->frag($frag);

  #error checking
  throw("Can't filter fastq without a run id") unless($self->run_id);
  throw("Max percent N can't be greater than 1") if($self->max_percent_n > 1);
  ######

  return $self;
}


=head2 get_input_files_from_collection

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : ReseqTrack::Collection
  Function  : To retrive the file paths from the collection object
  Returntype: arrayref of strings (filepaths)
  Exceptions: 
  Example   : 

=cut


sub get_input_files_from_collection{
  my ($self, $collection) = @_;
  $collection = $self->collection unless($collection);
  my $file_objects = $collection->others;
  my @array;
  foreach my $file_object(@$file_objects){
    push(@array, $file_object->name);
  }
  return \@array;
}

=head2 filter_files

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : arrayref, filepaths
  Function  : sort out which files are which and then filter them on the basis
  of the specified criteria
  Returntype: arrayref of compressed output files
  Exceptions: throws if files don't exist, throws if files were assigned to mate1, 
  mate2 or frag appropriately, throws if output files already exist and clobber not
  set
  Example   : 

=cut



sub filter_files{
  my ($self, $file_paths) = @_;
  #assign files to mate1, 2 and frag
  my $file_count;
  $file_count++ if($self->mate1);
  $file_count++ if($self->mate2);
  $file_count++ if($self->frag);
  unless($file_count){
    if(!$file_paths || @$file_paths == 0){
      $file_paths = $self->input_files;
      if(!$file_paths || @$file_paths == 0){
	$file_paths = $self->get_input_files_from_collection;
      }
    }
    my ($mate1, $mate2, $frag) = assign_files($file_paths);
    $self->mate1($mate1);
    $self->mate2($mate2);
    $self->frag($frag);
  }
  #throw if incorrect set of files present
  if(($self->mate1 && !(-e $self->mate1)) ||($self->mate2 && !(-e $self->mate2)) || ($self->frag && !(-e $self->frag))){
    throw("Don't have any filepaths to filter on");
  }
  if(($self->mate1 && !$self->mate2) || (!$self->mate1 && $self->mate2)){
    throw("Dont have both mate1 and mate2 from ".$self->run_id);
  }
  #if mate1 and 2 are defined, setup output paths and filter, ensuring the frags
  #go into a new frag file
  if($self->mate1 && $self->mate2){
    my $filtered_mate1_name = $self->create_filtered_filename(basename($self->mate1));
    my $path = $self->output_dir."/".$filtered_mate1_name;
    my $gz_path = $path.".gz";
    if(-s $path || -s $gz_path){
      if($self->clobber){
	delete_file($path) if(-e $path);
	delete_file($gz_path) if(-e $gz_path);
      }else{
	throw("Can't run with ".$path." already exists");
      }
    }
    $self->filtered_mate1($self->output_dir."/".$filtered_mate1_name);
    my $filtered_mate2_name = $self->create_filtered_filename(basename($self->mate2));
    $path = $self->output_dir."/".$filtered_mate1_name;
    $gz_path = $path.".gz";
    if(-s $path || -s $gz_path){
      if($self->clobber){
	delete_file($path) if(-e $path);
	delete_file($gz_path) if(-e $gz_path);
      }else{
	throw("Can't run with ".$path." already exists");
      }
    }
    $self->filtered_mate2($self->output_dir."/".$filtered_mate2_name);
  }
  #sorting out frag filtering
  my $filtered_frag_name = $self->create_filtered_filename(basename($self->frag)) 
    if($self->frag);
  if(!$filtered_frag_name){
    $filtered_frag_name = $self->run_id.".filt.fastq";
  }
  my $path = $self->output_dir."/".$filtered_frag_name;
  my $gz_path = $path.".gz";
  if(-s $path || -s $gz_path){
    if($self->clobber){
      delete_file($path) if(-e $path);
      delete_file($gz_path) if(-e $gz_path)
    }else{
      throw("Can't run with ".$path." already exists");
    }
  }
  #Filtering all files
  $self->filtered_frag($self->output_dir."/".$filtered_frag_name);
  $self->filter_matepairs if($self->mate1 && $self->mate2);
  $self->filter_frag if($self->frag);
  $self->close_files;
  #compressing files
  my @compressed_files;
  if($self->filtered_mate1 && -e $self->filtered_mate1){
    unless(-s $self->filtered_mate1){
      #print STDERR $self->filtered_mate1." is empty deleting\n";
      delete_file($self->filtered_mate1);
    }else{
      push(@compressed_files, $self->compress_file($self->filtered_mate1));
    }
  }
  if($self->filtered_mate2 && -e $self->filtered_mate2){
    unless(-s $self->filtered_mate2){
      #print STDERR $self->filtered_mate2." is empty deleting\n";
      delete_file($self->filtered_mate2);
    }else{
      push(@compressed_files, $self->compress_file($self->filtered_mate2));
    }
  }
  if($self->filtered_frag && -e $self->filtered_frag){
    unless(-s $self->filtered_frag){
      #print STDERR $self->filtered_frag." is empty deleting\n";
      delete_file($self->filtered_frag);
    }else{
      push(@compressed_files, $self->compress_file($self->filtered_frag));
    }
  }
  $self->compressed_files(\@compressed_files);
  $self->print_error_hash_summary();
  return \@compressed_files;
}




=head2 compress_file

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, filepath
  Function  : run gzip to compress the given file
  Returntype: string, path to gzipped file
  Exceptions: throw if non zero exit code or if fails to create the gzipped file
  Example   : 

=cut



sub compress_file{
  my ($self, $file) = @_;
  my $cmd = "gzip ".$file;
  my $new_file = $file.".gz";
  if(-e $new_file){
    unlink $new_file;
  }
  my $exit = system($cmd);
  if($exit >= 1){
    throw("There is a problem with ".$cmd." non zero exit code ".$exit);
  }
  unless(-e $new_file){
    throw("Failed to produce ".$new_file." from ".$file." using ".$cmd);
  }
  return $new_file;
}


=head2 filter_frag

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Function  : filter fragment file
  Returntype: 
  Exceptions: throws if string is colorspace but boolean is_colorspace is not set
  Example   : 

=cut



sub filter_frag{
  my ($self) = @_;
  #get file handles
  my $frag_fh = $self->frag_fh;
  my $filt_frag_fh = $self->filtered_frag_fh;
  my $unfiltered_readcount;
  my $unfiltered_basecount;
  #establish starting base and read counts
  my $frag_readcount = $self->filtered_frag_readcount;
  my $frag_basecount = $self->filtered_frag_basecount;
  my $run_id = $self->run_id;
  #start reading file
  while ( my $header = <$frag_fh> ) {
    my $seq = <$frag_fh>;
    my $strand = <$frag_fh>;
    my $qual = <$frag_fh>;
    if($self->looks_like_colorspace($seq) && !($self->is_colorspace)){
      throw("Have colorspace ".$header." ".$seq." but it isn't set to colorspace")
    }
    #remove new line characters
    chomp ($header, $seq, $strand, $qual);
    $header =~ s/^\s+|\s+$//g;
    $seq =~ s/^\s+|\s+$//g;
    $strand =~ s/^\s+|\s+$//g;
    $qual =~ s/^\s+|\s+$//g;
    #increment read and base count
    $unfiltered_basecount += $self->read_length($seq);
    $unfiltered_readcount++;
    #check the read
    my $status = $self->check_read($header, $seq, $strand, $qual, 
				   $self->frag, $self->run_id,
				   $unfiltered_readcount);
    if($status){
      #print read and alter base and read count
      $self->print_read($header, $seq, $strand, $qual, $filt_frag_fh);
      $frag_basecount += $self->read_length($seq);
      $frag_readcount++;
    }
  }
  #Update base and read counts
  $self->unfiltered_frag_basecount($unfiltered_basecount);
  $self->unfiltered_frag_readcount($unfiltered_readcount);
  $self->filtered_frag_readcount($frag_readcount);
  $self->filtered_frag_basecount($frag_basecount);
}


=head2 filter_matepairs

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Function  : iterate through mate1 and mate2 files ensuring all reads are of
  appropriate quality and they match in length etc
  Returntype: n/a
  Exceptions: throw if have colorspace read but is_colorspace is false
  Example   : 

=cut



sub filter_matepairs{
  my ($self) = @_;
  #get file handles
  my $mate1_fh = $self->mate1_fh;
  my $mate2_fh = $self->mate2_fh;
  my $filt_mate1_fh = $self->filtered_mate1_fh;
  my $filt_mate2_fh = $self->filtered_mate2_fh;
  my $filt_frag_fh = $self->filtered_frag_fh;
  #establish starting base and read counts
  my $unfiltered_basecount_1 = 0;
  my $unfiltered_basecount_2 = 0;
  my $filtered_basecount_1 = 0;
  my $filtered_basecount_2 = 0;
  my $frag_basecount = 0;
  my $unfiltered_readcount_1 = 0;
  my $unfiltered_readcount_2 = 0;
  my $filtered_readcount_1 = 0;
  my $filtered_readcount_2 = 0;
  my $frag_readcount = 0;
  my $run_id = $self->run_id;
  #start reading file
  while ( my $header1 = <$mate1_fh> ) {
    my $seq1 = <$mate1_fh>;
    my $strand1 = <$mate1_fh>;
    my $qual1 = <$mate1_fh>;
    if($self->looks_like_colorspace($seq1) && !($self->is_colorspace)){
      throw("Have colorspace ".$header1." ".$seq1." but it isn't set to colorspace")
    }
    my $header2 = <$mate2_fh>;
    my $seq2 = <$mate2_fh>;
    my $strand2 = <$mate2_fh>;
    my $qual2 = <$mate2_fh>;
    #remove new line characters
    chomp ($header1, $header2, $seq1, $seq2, $strand1, $strand2, $qual1, $qual2);
    $header1 =~ s/^\s+|\s+$//g;
    $seq1 =~ s/^\s+|\s+$//g;
    $strand1 =~ s/^\s+|\s+$//g;
    $qual1 =~ s/^\s+|\s+$//g;
    
    $header2 =~ s/^\s+|\s+$//g;
    $seq2 =~ s/^\s+|\s+$//g;
    $strand2 =~ s/^\s+|\s+$//g;
    $qual2 =~ s/^\s+|\s+$//g;
    #increment read and base count
    $unfiltered_basecount_1 += $self->read_length($seq1);
    $unfiltered_readcount_1++;
    $unfiltered_basecount_2 += $self->read_length($seq2);
    $unfiltered_readcount_2++;
    my $seq1_status = 1;
    my $seq2_status = 1;
    #check both reads
    $seq1_status = $self->check_read($header1, $seq1, $strand1, $qual1, 
				     $self->mate1, $self->run_id,
				     $unfiltered_readcount_1);
    $seq2_status = $self->check_read($header2, $seq2, $strand2, $qual2, 
				     $self->mate2, $self->run_id,
				     $unfiltered_readcount_2);
    #if both are good print both and increment appropriate counts
    if($seq1_status && $seq2_status){
      $self->print_read($header1, $seq1, $strand1, $qual1, $filt_mate1_fh);
      $self->print_read($header2, $seq2, $strand2, $qual2, $filt_mate2_fh);
      $filtered_basecount_1 += $self->read_length($seq1);
      $filtered_readcount_1++;
      $filtered_basecount_2 += $self->read_length($seq2);
      $filtered_readcount_2++;
    }else{
      #else print the appropriate one and the frag count
      if($seq1_status){
	$self->print_read($header1, $seq1, $strand1, $qual1, $filt_frag_fh);
	$frag_basecount += $self->read_length($seq1);
	$frag_readcount++;
      }elsif($seq2_status){
	$self->print_read($header2, $seq2, $strand2, $qual2, $filt_frag_fh);
	$frag_basecount += $self->read_length($seq2);
	$frag_readcount++;
      }
    }
  }
  #Update base and read counts
  $self->unfiltered_mate1_basecount($unfiltered_basecount_1);
  $self->unfiltered_mate1_readcount($unfiltered_readcount_1);
  $self->unfiltered_mate2_basecount($unfiltered_basecount_2);
  $self->unfiltered_mate2_readcount($unfiltered_readcount_2);
  $self->filtered_mate1_basecount($filtered_basecount_1);
  $self->filtered_mate1_readcount($filtered_readcount_1);
  $self->filtered_mate2_basecount($filtered_basecount_2);
  $self->filtered_mate2_readcount($filtered_readcount_2);
  $self->filtered_frag_readcount($frag_readcount);
  $self->filtered_frag_basecount($frag_basecount);
}


=head2 print_read

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, fastq header
  Arg [3]   : string, sequence
  Arg [4]   : string, strand
  Arg [5]   : string, qual
  Arg [6]   : filehandle
  Function  : print fastq data to filehandle
  Returntype: 1
  Exceptions: 
  Example   : 

=cut



sub print_read{
   my ($self, $header, $seq, $strand, $qual, $fh) = @_;
   if(!$fh){
     throw("Have no filehandle to print to");
   }
   unless(fileno($fh)){
     throw("Can't print ".$header." to an undefined filehandle");
   }
   print $fh $header."\n";
   print $fh $seq."\n";
   print $fh $strand."\n";
   print $fh $qual."\n";
   return 1;
}


=head2 check_read

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, fastq header
  Arg [3]   : string, sequence
  Arg [4]   : string, strand
  Arg [5]   : string, qual
  Arg [6]   : string, filename
  Arg [7]   : string, run id
  Arg [8]   : int, $read_count
  Function  : run through all checks for a read set
  Returntype: boolean
  Exceptions: throws if syntax errors are found
  Example   : 

=cut



sub check_read{
  my ($self, $header, $seq, $strand, $qual, $filename, $run_id, $read_count) = @_;
  #print STDERR "Filtering".$header."\n";
  if(!$header || !$seq || !$strand || !$qual){
    #print STDERR "Header ".$header."\n";
    #print STDERR "SEQ ".$seq."\n";
    #print STDERR "strand ".$strand."\n";
    #print STDERR "qual ".$qual."\n";
    #throw($filename." ".$run_id." ".$read_count." has an empty string\n");
    #print "Empty strings returning 0\n";
    return 0;
  }

  # should not have non-printing characters in seq or qual lines
  return 0 if  ( $self->non_printing_char_check ($header, $seq , "seq_line") );
  return 0 if  ( $self->non_printing_char_check ($header, $qual, "qual_line") );


  $run_id = $self->run_id unless($run_id);
  unless($header =~ /^\@/){
    throw($header." from ".$filename." does not start with @ like it should");
  }
  unless($strand =~ /^\+/){
    throw($strand." from ".$filename." does not start with +  like it should");
  }
  unless($header =~ /$run_id/){
    #print STDERR $header." from ".$filename." does not contain the run id\n";
    #should die here
    $self->add_to_error_hash("Run id missing from header");
    return 0;
  }
  if(length($strand) > 1 && !($strand =~ /$run_id/)){
    #print STDERR $header." ".$strand." from ".$filename." does not contain ".
    #  "the run id\n";
    #should die here
    $self->add_to_error_hash("Run id missing from strand string");
    return 0;
  }
  unless(length($seq) == length($qual)){
    #print STDERR "Seq and Qual strings from ".$filename." ".$header.
    #  " are of different lengths\n";
    $self->add_to_error_hash("Sequence and Qual strings different lengths");
    return 0;
  }
  unless($self->check_length($seq)){
    #print STDERR "There is a problem with the length of ".$filename." ".$read_count.
    #  " seq or qual string ".$header."\n";
    $self->add_to_error_hash("sequence string too short");
    return 0;
  }
  unless($self->check_prop_n($seq)){
    #print STDERR "Sequence string  for ".$filename." ".$read_count.
    #  " contains too many Ns ".$header." \n"; 
    $self->add_to_error_hash("Too Many Ns");
    return 0;
  }
  #print "Checking ".$header."\n";
  #print $seq."\n";
  unless($self->check_base_comp($seq, $header)){
    #print STDERR $header."\n";
    #print STDERR "Sequence string  for ".$filename." ".$read_count.
    #  " contains too many runs of the same base ".$header."\n"; 
    $self->add_to_error_hash("Too many bases of same type");
    return 0;
  }
  unless($self->check_qual($qual)){
    #print STDERR "Qual string  for ".$filename." ".$read_count.
    #  " contains too many low values ".$header."\n"; 
    $self->add_to_error_hash("Too low qual values");
    return 0;
  }
  return 1;
}



=head2 check_prop_n

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, sequence string
  Function  : checks if there are too many Ns or .s in the sequence
  Returntype: boolean
  Exceptions: 
  Example   : 

=cut



sub check_prop_n{
  my ($self, $seq) = @_;
  my @chars = split //, $seq;
  my $length = length($seq);
  my $count = 0;
  if($self->is_colorspace){
    $length -= 1;
    shift @chars;
    foreach my $char(@chars){
      $count++ if($char eq '.' || $char eq 'N');
    }
  }else{
    foreach my $char(@chars){
      $count++ if($char eq 'N');
    }
  }

  if(($count/$length) > $self->max_percent_n){
    #print STDERR "Have more than 50% n \n";
    return 0;
  }
  return 1;
}



=head2 looks_like_colorspace

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, sequence string
  Function  : checks if the string it colorspace or not
  Returntype: boolean
  Exceptions: 
  Example   : 

=cut



sub looks_like_colorspace{
  my ($self, $string) = @_;
  return 0 if($string =~ /[ATGCN]/);
  return 1 if($string =~ /[0123\.N]/);
  return 0;
}


=head2 check_base_comp

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, sequence string
  Function  : check that the sequence string isn't to many Ns or another single
  base pair
  Returntype: boolean 
  Exceptions: 
  Example   : 

=cut



sub check_base_comp{
  my ($self, $seq, $header) = @_;
  my $len = $self->min_length;
  my ($seq_str) = $self->get_seq_string_to_check($seq);
  if($seq =~ /([^A|C|G|T|N|0|1|2|3|\.|\n])/gi){
    #print "Throwing away ".$header." ".$seq."\n";
    #throw("Non standard characters in ".$header." ".$seq);
    return 0;
  }
  if($self->is_colorspace){
    my $length = length($seq_str) - 1;
    my $colour_space_string = substr($seq, 1, $length); 
    if($colour_space_string =~ /^[A-Z]/i
	    || !($self->looks_like_colorspace($colour_space_string))) {
      throw("There is a problem for ".$header." ".$seq." colorspace run as ".
	    "still contains letters or non numeric characters after subseq ".
	    $seq_str);
      return 0;
    }
    if($seq_str =~ /\./ || $seq_str =~ /N/i){
      #print $seq_str." is wrong as it contains . or N\n";
      return 0;
    }elsif($seq_str =~ /0{$len}/i){
      return 0;
    }elsif($seq_str =~ /1{$len}/i){
      return 0;
    }elsif($seq_str =~ /2{$len}/i){
      return 0;
    }elsif($seq_str =~ /3{$len}/i){
      return 0;
    }
    #print "Colorspace returning true\n";
    return 1;
  }else{
    if ($seq_str =~ /N/i) {
      #print  "ERROR: in Block $seq_str, the sequence contains Ns in the first $len bp region\n";
      return 0;
    }
    elsif ($seq_str =~ /A{$len}/i ){ 
     # print  "ERROR: in Block $seq_str, the sequence contains only AAAAAAAAAAAAA in the first $len bp\n";
      return 0;
    }
    elsif ($seq_str =~ /T{$len}/i ) {
      #print  "ERROR: in Block $seq_str, the sequence contains only TTTTTTTTTTTT in the first $len bp\n";
      return 0;
    }
    elsif ($seq_str =~ /G{$len}/i ) {
      #print  "ERROR: in Block $seq_str, the sequence contains only GGGGGGGGGGGGGG in the first $len bp\n";
      return 0;
    }
    elsif ($seq_str =~ /C{$len}/i ) {
      #print  "ERROR: in Block $seq_str, the sequence contains only CCCCCCCCCCCCCCC in the first $len bp\n";
      return 0;
    }
    return 1;
  }
}


=head2 check_length

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, seq string
  Function  : check if sequence string is longer than minimum length
  Returntype: boolean
  Exceptions: 
  Example   : 

=cut


sub check_length{
  my ($self, $seq) = @_;
  if(!$seq){
    throw("Have no seq string to check the length of");
  }
  return 0 if($self->read_length($seq) < $self->min_length);
  return 1;
}


=head2 check_qual

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, qual string
  Function  : checks that there aren't quality values less than minimum in the 
  first so many bases in the qual string
  Returntype: boolean
  Exceptions: 
  Example   : 

=cut



sub check_qual{
  my ($self, $qual) = @_;
  my ($qual_str) = $self->get_qual_string_to_check($qual);
  my @values = split //, $qual_str;
  my $count = 0;
  foreach my $value(@values){
    my $q_num = ord($value) - 33;
    #print STDERR "Checking ".$value." ".$q_num."\n";
    $count++ if($q_num < $self->min_qual);
  }
  return 0 if($count);
  return 1;
}


=head2 check_string_contents

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, sequence string
  Function  : checks if sequence string has any unexpected characters in it
  Returntype: boolean
  Exceptions: 
  Example   : 

=cut



sub check_string_contents{
  my ($self, $seq) = @_;
  my $substr;
  if($self->is_colorspace){
    $substr = substr($seq, 1);
    return $self->looks_like_colorspace($substr);
  }else{
    return $self->check_base_seq($seq);
  }
}

=head2 check_base_seq

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, seq string
  Function  : ensures the sequence is only made up of bases ATGCN
  Returntype: boolean
  Exceptions: 
  Example   : 

=cut


sub check_base_seq{
  my ($self, $seq) = @_;
  return 0 unless($seq =~ /^[ATGCN]$/i);
  return 1;
}

=head2 get_qual_string_to_check

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, qual string
  Function  : returns the minimum number of bases from the qual string
  Returntype: string
  Exceptions: 
  Example   : 

=cut


sub get_qual_string_to_check{
  my ($self, $qual) = @_;
  if(!$qual){
    throw("Have no seq or qual string to check");
  }
  my $qual_offset = 0;
  my $qual_length = $self->min_length;
  if($self->is_colorspace){
    $qual_offset = 1;
  }
  my $qual_substr = substr($qual, $qual_offset, $qual_length);
  return $qual_substr;
}


=head2 get_seq_string_to_check

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, seq string
  Function  : returns the minimum number of bases from the seq string
  Returntype: string
  Exceptions: 
  Example   : 

=cut



sub get_seq_string_to_check{
  my ($self, $seq) = @_;
  if(!$seq){
    throw("Have no seq string to check");
  }
  my $seq_offset = 0;
  my $seq_length = $self->min_length;
  if($self->is_colorspace){
    $seq_offset = 1;
    $seq_length = $seq_length+1;
  }
  my $seq_substr = substr($seq, $seq_offset, $seq_length);
  return ($seq_substr);
}


=head2 read_length

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, sequence string
  Function  : calculates length of string
  Returntype: int
  Exceptions: n/a
  Example   : my $length = $self->read_length($seq);

=cut



sub read_length{
  my ($self, $seq) = @_;
  if($self->is_colorspace){
    return length($seq) -1;
  }else{
    return length($seq);
  }
}


=head2 close_files

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Function  : undef all possibly open file handles
  Returntype: n/a
  Exceptions: n/a
  Example   : $self->close_files();

=cut



sub close_files{
  my ($self) = @_;
  $self->{mate1_fh} = undef;
  $self->{mate2_fh} = undef;
  $self->{frag_fh} = undef;
  $self->{filtered_mate1_fh} = undef;
  $self->{filtered_mate2_fh} = undef;
  $self->{filtered_frag_fh} = undef;
}



=head2 create_filtered_filename

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, filename
  Function  : convert a archive fastq name (nme.fastq.gz) to a filtered fastq name
  (name.filt.fastq);
  Returntype: string
  Exceptions: n/a
  Example   : my $filtered_mate1 = $self->create_filtered_filename($mate1);

=cut



sub create_filtered_filename{
  my ($self, $name) = @_;
  $name =~ s/\.fastq\.gz/\.filt\.fastq/;
  return $name;
}



=head2 create_fh

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string, filepath
  Arg [3]   : string, mode (program to use to open file) or read write mode
  Arg [4]   : boolean, flag to indicate | is needed at end of open command
  Function  : to open a filehandle for the process, it can use a non standard
  open command if told, by default it will use zcat to open gz files unless otherwise
  specified
  Returntype: FileHandle 
  Exceptions: throws if it failed to open file
  Example   : my $fh = $self->create_fh($self->mate1);

=cut



sub create_fh{
  my ($self, $file, $mode, $use_pipe) = @_;
  my $fh = new FileHandle;
  throw("Can't process an empty string as a file") unless($file);
  if($file =~ /\.gz$/){
    $mode = "zcat" unless($mode);
    $use_pipe = 1 unless(defined($use_pipe));
  }
  my $open_cmd = $mode." ".$file;
  $open_cmd .= " | " if($use_pipe);
  eval{
    $fh->open($open_cmd);
  };
  if($@){
    throw("Problem with ".$open_cmd." $@");
  }
  if(!$fh || !fileno($fh)){
    throw("Failed to open ".$open_cmd);
  }
  return $fh;
}


sub print_error_hash_summary{
  my ($self) = @_;
  print STDERR "\n";
  foreach my $key(keys(%{$self->error_hash})){
    print STDERR $key." ".$self->{error_hash}->{$key}."\n";
  }
  print STDERR "\n";
}
sub add_to_error_hash{
  my ($self, $key) = @_;
  throw("Can't add an empty key to hash") unless($key);
  unless($self->{error_hash}->{$key}){
    $self->{error_hash}->{$key} = 1;
  }else{
    $self->{error_hash}->{$key}++;
  }
}

sub error_hash{
  my ($self, $arg) = @_;
  if($arg){
    $self->{error_hash} = $arg;
  }
  if(!$self->{error_hash}){
    $self->{error_hash} = {};
  }
  return $self->{error_hash};
}


=head2 non_printing_char_check
  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Arg [2]   : string to check for nonprinting chars
  Arg [3]   : header of read
  Arg [4]   : line type ( seq or qual)
  Function  : check for characters in ascii range 3-127
  Returntype: boolean
  Exceptions: NULL string passed
  Example   :   $self->non_printing_char_check ($seq ,$header, "seq_line");

=cut

sub  non_printing_char_check{

  my ( $self, $header, $input_string , $descriptor)= @_;

  throw ("NULL string passed\n") if (!defined $input_string);

  $input_string =~ /([\x80-\xFF]|[\x00-\x20])/;

  if ( $1 ){
    my $bad_char_ord = ord($1);
    my $msg = "read: $header ";
    $msg .= "contains nonprinting character $descriptor (ascii = $bad_char_ord)\n";
    print $msg;
    $self->add_to_error_hash($msg);
    return 1;
  }

 return 0;
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

sub compressed_files{
  my ($self, $arg) = @_;
  if($arg){
    $self->{compressed_files} = $arg;
  }
  return $self->{compressed_files};
}

sub output_dir{
  my ($self, $output_dir) = @_;
  if($output_dir){
    throw($output_dir." must exist as a directory") unless(-d $output_dir);
    $self->{output_dir} = $output_dir;
  }
  $self->{output_dir};
}

sub collection{
  my ($self, $collection) = @_;
  if($collection){
    assert_ref($collection, 'ReseqTrack::Collection');
    $self->{collection} = $collection;
  }
  return $self->{collection};
}

sub min_length{
  my ($self, $min_length) = @_;
  if($min_length){
    $self->{min_length} = $min_length;
  }
  return $self->{min_length};
}
sub min_qual{
  my ($self, $min_qual) = @_;
  if($min_qual){
    $self->{min_qual} = $min_qual;
  }
  return $self->{min_qual};
}

sub max_percent_n{
  my ($self, $max_percent_n) = @_;
  if($max_percent_n){
    $self->{max_percent_n} = $max_percent_n;
  }
  return $self->{max_percent_n};
}

sub is_colorspace{
  my ($self, $is_colorspace) = @_;
  if($is_colorspace){
    $self->{is_colorspace} = $is_colorspace;
  }
  return $self->{is_colorspace};
}

sub run_id{
  my ($self, $run_id) = @_;
  if($run_id){
    $self->{run_id} = $run_id;
  }
  return $self->{run_id};
}

sub unfiltered_mate1_basecount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{unfiltered_mate1_basecount} = $count;
  }
  return $self->{unfiltered_mate1_basecount};
}

sub unfiltered_mate1_readcount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{unfiltered_mate1_readcount} = $count;
  }
  return $self->{unfiltered_mate1_readcount};
}

sub unfiltered_mate2_basecount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{unfiltered_mate2_basecount} = $count;
  }
  return $self->{unfiltered_mate2_basecount};
}

sub unfiltered_mate2_readcount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{unfiltered_mate2_readcount} = $count;
  }
  return $self->{unfiltered_mate2_readcount};
}

sub filtered_mate1_basecount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{filtered_mate1_basecount} = $count;
  }
  return $self->{filtered_mate1_basecount};
}

sub filtered_mate1_readcount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{filtered_mate1_readcount} = $count;
  }
  return $self->{filtered_mate1_readcount};
}

sub unfiltered_frag_basecount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{unfiltered_frag2_basecount} = $count;
  }
  return $self->{unfiltered_frag2_basecount};
}

sub unfiltered_frag_readcount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{unfiltered_frag2_readcount} = $count;
  }
  return $self->{unfiltered_frag2_readcount};
}

sub filtered_frag_basecount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{filtered_frag1_basecount} = $count;
  }
  return $self->{filtered_frag1_basecount};
}

sub filtered_frag_readcount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{filtered_frag1_readcount} = $count;
  }
  return $self->{filtered_frag1_readcount};
}


sub filtered_mate2_basecount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{filtered_mate2_basecount} = $count;
  }
  return $self->{filtered_mate2_basecount};
}

sub filtered_mate2_readcount{
  my ($self, $count) = @_;
  if(defined($count)){
    $self->{filtered_mate2_readcount} = $count;
  }
  return $self->{filtered_mate2_readcount};
}


sub mate1{
  my ($self, $mate1) = @_;
  if($mate1){
    throw("Can't run with a file that doesn't exist ".$mate1) unless(-e $mate1);
    $self->{mate1} = $mate1;
  }
  return $self->{mate1};
}

sub mate1_fh{
  my ($self, $fh) = @_;
  if($fh){
    $self->{mate1_fh} = $fh;
  }
  unless($self->{mate1_fh}){
    $fh = $self->create_fh($self->mate1);
    $self->{mate1_fh} = $fh;
  }
  return $self->{mate1_fh};
}

sub mate2{
  my ($self, $mate2) = @_;
  if($mate2){
    throw("Can't run with a file that doesn't exist ".$mate2) unless(-e $mate2);
    $self->{mate2} = $mate2;
  }
  return $self->{mate2};
}

sub mate2_fh{
  my ($self, $fh) = @_;
  if($fh){
    $self->{mate2_fh} = $fh;
  }
  unless($self->{mate2_fh}){
    $fh = $self->create_fh($self->mate2);
    $self->{mate2_fh} = $fh;
  }
  return $self->{mate2_fh};
}


sub frag{
  my ($self, $frag) = @_;
  if($frag){
    throw("Can't run with a file that doesn't exist ".$frag) unless(-e $frag);
    $self->{frag} = $frag;
  }
  return $self->{frag};
}


sub frag_fh{
  my ($self, $fh) = @_;
  if($fh){
    $self->{frag_fh} = $fh;
  }
  unless($self->{frag_fh}){
    $fh = $self->create_fh($self->frag);
    $self->{frag_fh} = $fh;
  }
  return $self->{frag_fh};
}

sub input_files{
  my ($self, $arrayref) = @_;
  if($arrayref){
    throw("Failed to give FilterFastq:input_files an arrayref but ".$arrayref)
      unless(ref($arrayref) eq 'ARRAY');
    $self->{input_files} = $arrayref;
  }
  return $self->{input_files};
}

sub filtered_mate1{
  my ($self, $filtered_mate1) = @_;
  if($filtered_mate1){
     $self->{filtered_mate1} = $filtered_mate1;
  }
  return $self->{filtered_mate1};
}

sub filtered_mate1_fh{
  my ($self, $fh) = @_;
  if($fh){
    $self->{filtered_mate1_fh} = $fh;
  }
  unless($self->{filtered_mate1_fh}){
    #print STDERR "Opening filtered mate1 file handle ".$self->filtered_mate1."\n";
    $fh = $self->create_fh($self->filtered_mate1, ">") if($self->filtered_mate1);
    $self->{filtered_mate1_fh} = $fh;
  }
  return  $self->{filtered_mate1_fh};
}

sub filtered_mate2{
  my ($self, $filtered_mate2) = @_;
  if($filtered_mate2){
     $self->{filtered_mate2} = $filtered_mate2;
  }
  return $self->{filtered_mate2};
}

sub filtered_mate2_fh{
  my ($self, $fh) = @_;
  if($fh){
    $self->{filtered_mate2_fh} = $fh;
  }
  unless($self->{filtered_mate2_fh}){
    $fh = $self->create_fh($self->filtered_mate2, ">") if($self->filtered_mate2);
    $self->{filtered_mate2_fh} = $fh;
  }
  return  $self->{filtered_mate2_fh};
}


sub filtered_frag{
  my ($self, $filtered_frag) = @_;
  if($filtered_frag){
     $self->{filtered_frag} = $filtered_frag;
  }
  return $self->{filtered_frag};
}

sub filtered_frag_fh{
  my ($self, $fh) = @_;
  if($fh){
    $self->{filtered_frag_fh} = $fh;
  }
  unless($self->{filtered_frag_fh}){
    $fh = $self->create_fh($self->filtered_frag, ">>") if($self->filtered_frag);
    $self->{filtered_frag_fh} = $fh;
  }
  return $self->{filtered_frag_fh};
}

sub clobber{
  my ($self, $clobber) = @_;
  if(defined($clobber)){
    $self->{clobber} = $clobber;
  }
  return $self->{clobber};
}

1;
