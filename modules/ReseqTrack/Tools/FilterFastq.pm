=pod

=head1 NAME

ReseqTrack::Tools::FilterFastq

=head1 SYNOPSIS

This is a object to filter fastq files on the basis of length of read, quality values
and percent Ns

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
use FileHandle;
use File::Basename;

=head2 new

  Arg [1]   : ReseqTrack::Tools::FilterFastq
  Function  : 
  Returntype: ReseqTrack::Tools::FilterFastq
  Exceptions: 
  Example   : 

=cut


sub new {
  my ($class, @args) = @_;
  my $self ={};
  bless $self,$class;
  my ($collection, $rmi, $min_length, $min_qual, $max_percent_n, $output_dir, 
      $is_colorspace, $run_id) 
    = rearrange([qw(COLLECTION RMI MIN_LENGTH MIN_QUAL MAX_PERCENT_N
		    OUTPUT_DIR IS_COLORSPACE RUN_ID)],  @args);

  #SETTING DEFAULTS
  $self->min_length(25);
  $self->min_qual(2);
  $self->max_percent_n(50);
  $run_id = $rmi->run_id if(!$run_id && $rmi);
  ########

  $self->collection($collection);
  $self->rmi($rmi);
  $self->min_length($min_length);
  $self->min_qual($min_qual);
  $self->max_percent_n($max_percent_n);
  $self->is_colorspace($is_colorspace);
  $self->output_dir($output_dir);
  $self->run_id($run_id);

  #error checking
  throw("Can't filter fastq files without a collection object") unless($self->collection);
  throw("Can't filter fastq without a run id") unless($self->run_id);
  ######

  return $self;
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

sub rmi{
  my ($self, $rmi) = @_;
  if($rmi){
    assert_ref($rmi, 'ReseqTrack::RunMetaInfo');
    $self->{rmi} = $rmi;
  }
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

sub filter_files{
  my ($self, $file_paths) = @_;
  unless($file_paths && @$file_paths >= 1){
    my $others = $self->collection->others;
    foreach my $file(@$others){
      push(@$file_paths, $file->name);
    }
  }
  my ($mate1, $mate2, $frag) = assign_files($file_paths);
  if(($mate1 && !(-e $mate1)) ||($mate2 && !(-e $mate2)) || ($frag && !(-e $frag))){
    print STDERR "Looking at ".$self->collection->name."\n";
    throw("Don't have any filepaths to filter on");
  }
  if(($mate1 && !$mate2) || (!$mate1 && $mate2)){
    throw("Dont have both mate1 and mate2 from ".$self->collection->name);
  }
  #need writeable file handles for mate1 mate2 and fragment output
  if($mate1 && $mate2){
    $self->mate1($mate1);
    my $filtered_mate1_name = $self->create_filtered_filename(basename($mate1));
    $self->filtered_mate1($self->output_dir."/".$filtered_mate1_name);
    $self->mate2($mate2);
    my $filtered_mate2_name = $self->create_filtered_filename(basename($mate2));
    $self->filtered_mate2($self->output_dir."/".$filtered_mate2_name);
  }
  $self->frag($frag);
  my $filtered_frag_name = $self->create_filtered_filename(basename($frag)) if($frag);
  if(!$filtered_frag_name){
    $filtered_frag_name = $self->collection->name.".filt.fastq";
  }
  $self->filtered_frag($self->output_dir."/".$filtered_frag_name);
  $self->filter_matepairs if($mate1 && $mate2);
  $self->filter_frag if($frag);
  $self->close_files;
  my @compressed_files;
  if($self->filtered_mate1 && -e $self->filtered_mate1){
    push(@compressed_files, $self->compress_file($self->filtered_mate1));
  }
  if($self->filtered_mate2 && -e $self->filtered_mate2){
    push(@compressed_files, $self->compress_file($self->filtered_mate2));
  }
  if($self->filtered_frag && -e $self->filtered_frag){
    push(@compressed_files, $self->compress_file($self->filtered_frag));
  }
  $self->compressed_files(\@compressed_files);
  return \@compressed_files;
}

sub compressed_files{
  my ($self, $arg) = @_;
  if($arg){
    $self->{compressed_files} = $arg;
  }
  return $self->{compressed_files};
}

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

sub filter_frag{
  my ($self) = @_;
  my $frag_fh = $self->frag_fh;
  my $filt_frag_fh = $self->filtered_frag_fh;
  my $unfiltered_readcount;
  my $unfiltered_basecount;
  my $frag_readcount = $self->filtered_frag_readcount;
  my $frag_basecount = $self->filtered_frag_basecount;
  my $run_id = $self->run_id;
  while ( my $header = <$frag_fh> ) {
    my $seq = <$frag_fh>;
    my $strand = <$frag_fh>;
    my $qual = <$frag_fh>;
    chomp ($header, $seq, $strand, $qual);
    $header =~ s/^\s+|\s+$//g;
    $seq =~ s/^\s+|\s+$//g;
    $strand =~ s/^\s+|\s+$//g;
    $qual =~ s/^\s+|\s+$//g;
    $unfiltered_basecount += $self->read_length($seq);
    $unfiltered_readcount++;
    my $status = $self->check_read($header, $seq, $strand, $qual, 
				   $self->frag, $self->run_id,
				   $unfiltered_readcount);
    if($status){
      $self->print_read($header, $seq, $strand, $qual, $filt_frag_fh);
      $frag_basecount += $self->read_length($seq);
      $frag_readcount++;
    }
  }
  $self->unfiltered_frag_basecount($unfiltered_basecount);
  $self->unfiltered_frag_readcount($unfiltered_readcount);
  $self->filtered_frag_readcount($frag_readcount);
  $self->filtered_frag_basecount($frag_basecount);
}

sub filter_matepairs{
  my ($self) = @_;
  my $mate1_fh = $self->mate1_fh;
  my $mate2_fh = $self->mate2_fh;
  my $filt_mate1_fh = $self->filtered_mate1_fh;
  my $filt_mate2_fh = $self->filtered_mate2_fh;
  my $filt_frag_fh = $self->filtered_frag_fh;
  my $unfiltered_basecount_1;
  my $unfiltered_basecount_2;
  my $filtered_basecount_1;
  my $filtered_basecount_2;
  my $frag_basecount;
  my $unfiltered_readcount_1;
  my $unfiltered_readcount_2;
  my $filtered_readcount_1;
  my $filtered_readcount_2;
  my $frag_readcount;
  my $run_id = $self->run_id;
  while ( my $header1 = <$mate1_fh> ) {
    my $seq1 = <$mate1_fh>;
    my $strand1 = <$mate1_fh>;
    my $qual1 = <$mate1_fh>;
    
    my $header2 = <$mate2_fh>;
    my $seq2 = <$mate2_fh>;
    my $strand2 = <$mate2_fh>;
    my $qual2 = <$mate2_fh>;
    chomp ($header1, $header2, $seq1, $seq2, $strand1, $strand2, $qual1, $qual2);
    $header1 =~ s/^\s+|\s+$//g;
    $seq1 =~ s/^\s+|\s+$//g;
    $strand1 =~ s/^\s+|\s+$//g;
    $qual1 =~ s/^\s+|\s+$//g;
    
    $header2 =~ s/^\s+|\s+$//g;
    $seq2 =~ s/^\s+|\s+$//g;
    $strand2 =~ s/^\s+|\s+$//g;
    $qual2 =~ s/^\s+|\s+$//g;
    $unfiltered_basecount_1 += $self->read_length($seq1);
    $unfiltered_readcount_1++;
    $unfiltered_basecount_2 += $self->read_length($seq2);
    $unfiltered_readcount_2++;
    my $seq1_status = 1;
    my $seq2_status = 1;
    $seq1_status = $self->check_read($header1, $seq1, $strand1, $qual1, 
				     $self->mate1, $self->run_id,
				     $unfiltered_readcount_1);
    $seq2_status = $self->check_read($header2, $seq2, $strand2, $qual2, 
				     $self->mate2, $self->run_id,
				     $unfiltered_readcount_2);
    #print $header1." SEQ1 status ".$seq1_status." SEQ2 status ".$seq2_status."\n";
    if($seq1_status && $seq2_status){
      #print "PRINTING BOTH\n";
      $self->print_read($header1, $seq1, $strand1, $qual1, $filt_mate1_fh);
      $self->print_read($header2, $seq2, $strand2, $qual2, $filt_mate2_fh);
      $filtered_basecount_1 += $self->read_length($seq1);
      $filtered_readcount_1++;
      $filtered_basecount_2 += $self->read_length($seq2);
      $filtered_readcount_2++;
    }else{
      if($seq1_status){
	$self->print_read($header1, $seq1, $strand1, $qual1, $filt_frag_fh);
	$frag_basecount += $self->read_length($seq1);
	$frag_readcount++;
      }elsif($seq2_status){
	$self->print_read($header2, $seq2, $strand2, $qual2, $filt_frag_fh);
	$frag_basecount += $self->read_length($seq2);
	$frag_readcount++;
      }else{
	#print STDERR "Skipping both reads for ".$run_id." ".
	#  $unfiltered_readcount_1."\n";
      }
    }
  }
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


sub unfiltered_mate1_basecount{
  my ($self, $count) = @_;
  if($count){
    $self->{unfiltered_mate1_basecount} = $count;
  }
  return $self->{unfiltered_mate1_basecount};
}

sub unfiltered_mate1_readcount{
  my ($self, $count) = @_;
  if($count){
    $self->{unfiltered_mate1_readcount} = $count;
  }
  return $self->{unfiltered_mate1_readcount};
}

sub unfiltered_mate2_basecount{
  my ($self, $count) = @_;
  if($count){
    $self->{unfiltered_mate2_basecount} = $count;
  }
  return $self->{unfiltered_mate2_basecount};
}

sub unfiltered_mate2_readcount{
  my ($self, $count) = @_;
  if($count){
    $self->{unfiltered_mate2_readcount} = $count;
  }
  return $self->{unfiltered_mate2_readcount};
}

sub filtered_mate1_basecount{
  my ($self, $count) = @_;
  if($count){
    $self->{filtered_mate1_basecount} = $count;
  }
  return $self->{filtered_mate1_basecount};
}

sub filtered_mate1_readcount{
  my ($self, $count) = @_;
  if($count){
    $self->{filtered_mate1_readcount} = $count;
  }
  return $self->{filtered_mate1_readcount};
}

sub unfiltered_frag_basecount{
  my ($self, $count) = @_;
  if($count){
    $self->{unfiltered_frag2_basecount} = $count;
  }
  return $self->{unfiltered_frag2_basecount};
}

sub unfiltered_frag_readcount{
  my ($self, $count) = @_;
  if($count){
    $self->{unfiltered_frag2_readcount} = $count;
  }
  return $self->{unfiltered_frag2_readcount};
}

sub filtered_frag_basecount{
  my ($self, $count) = @_;
  if($count){
    $self->{filtered_frag1_basecount} = $count;
  }
  return $self->{filtered_frag1_basecount};
}

sub filtered_frag_readcount{
  my ($self, $count) = @_;
  if($count){
    $self->{filtered_frag1_readcount} = $count;
  }
  return $self->{filtered_frag1_readcount};
}


sub filtered_mate2_basecount{
  my ($self, $count) = @_;
  if($count){
    $self->{filtered_mate2_basecount} = $count;
  }
  return $self->{filtered_mate2_basecount};
}

sub filtered_mate2_readcount{
  my ($self, $count) = @_;
  if($count){
    $self->{filtered_mate2_readcount} = $count;
  }
  return $self->{filtered_mate2_readcount};
}

sub print_read{
   my ($self, $header, $seq, $strand, $qual, $fh) = @_;
   print $fh $header."\n";
   print $fh $seq."\n";
   print $fh $strand."\n";
   print $fh $qual."\n";
   return 1;
}

sub check_read{
  my ($self, $header, $seq, $strand, $qual, $filename, $run_id, $read_count) = @_;
  $run_id = $self->run_id unless($run_id);
  unless($header =~ /^\@/){
    throw($header." from ".$filename." does not start with @ like it should");
  }
  unless($strand =~ /^\+/){
    throw($strand." from ".$filename." does not start with +  like it should");
  }
  unless($header =~ /$run_id/){
    #print STDERR $header." from ".$filename." does not contain the run id\n";
    return 0;
  }
  if(length($strand) > 1 && !($strand =~ /$run_id/)){
    #print STDERR $header." ".$strand." from ".$filename." does not contain ".
    #  "the run id\n";
    return 0;
  }
  unless(length($seq) == length($qual)){
    #print STDERR "Seq and Qual strings from ".$filename." ".$header.
    #  " are of different lengths\n";
    return 0;
  }
  unless($self->check_length($seq)){
    #print STDERR "There is a problem with the length of ".$filename." ".$read_count.
    #  " seq or qual string ".$header."\n";
    return 0;
  }
  unless($self->check_prop_n($seq)){
    #print STDERR "Sequence string  for ".$filename." ".$read_count.
    #  " contains to many Ns ".$header." \n"; 
    return 0;
  }
  unless($self->check_base_comp($seq)){
    #print STDERR "Sequence string  for ".$filename." ".$read_count.
    #  " contains to many runs of the same base ".$header."\n"; 
    return 0;
  }
  unless($self->check_qual($qual)){
    #print STDERR "Qual string  for ".$filename." ".$read_count.
    #  " contains to many low values ".$header."\n"; 
    return 0;
  }
  return 1;
}


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

  if(($count/$length) > 0.5){
    #print STDERR "Have more than 50% n \n";
    return 0;
  }
  return 1;
}


sub check_base_comp{
  my ($self, $seq) = @_;
  my $len = $self->min_length;
  my ($seq_str) = $self->get_seq_string_to_check($seq);
  if($self->is_colorspace){
    unless ($seq =~ /^[A-Z]/i && looks_like_number($seq_str) ) {
      return 0;
    }
    if($seq_str =~ /\./ || $seq_str =~ /N/i){
      return 0;
    }
    return 1;
  }else{
    if ($seq_str =~ /N/i) {
      #print  "ERROR: in Block $seq_str, the sequence contains Ns in the first $len bp region\n";
      return 0;
    }
    elsif ($seq_str =~ /A{$len}/i ){ 
      #print  "ERROR: in Block $seq_str, the sequence contains only AAAAAAAAAAAAA in the first $len bp\n";
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
sub check_length{
  my ($self, $seq) = @_;
  my ($seq_str) = $self->get_seq_string_to_check($seq);
  return 0 if(length($seq_str) < $self->min_length);
  return 1;
}

sub check_qual{
  my ($self, $qual) = @_;
  my ($qual_str) = $self->get_qual_string_to_check($qual);
  my @values = split //, $qual_str;
  my $count = 0;
  foreach my $value(@values){
    my $q_num = ord($value) - 33;
    $count++ if($q_num < 2);
  }
  return 0 if($count >= 2);
  return 1;
}
sub check_string_contents{
  my ($self, $seq) = @_;
  my $substr;
  if($self->is_colorspace){
    $substr = substr($seq, 1);
    return looks_like_number($substr);
  }else{
    return $self->check_base_seq($seq);
  }
}

sub check_base_seq{
  my ($self, $seq) = @_;
  return 0 unless($seq =~ /^[ATGCN]$/i);
  return 1;
}

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

sub get_seq_string_to_check{
  my ($self, $seq) = @_;
  if(!$seq){
    throw("Have no seq or qual string to check");
  }
  my $seq_offset = 0;
  my $seq_length = $self->min_length;
  if($self->is_colorspace){
    $seq_length = $seq_length+1;
  }
  my $seq_substr = substr($seq, $seq_offset, $seq_length);
  return ($seq_substr);
}

sub read_length{
  my ($self, $seq) = @_;
  if($self->is_colorspace){
    return length($seq) -1;
  }else{
    return length($seq);
  }
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
    my $fh = $self->create_fh($self->mate1, "zcat", 1);
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
    my $fh = $self->create_fh($self->mate2, "zcat", 1);
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
    my $fh = $self->create_fh($self->frag, "zcat", 1);
    $self->{frag_fh} = $fh;
  }
  return $self->{frag_fh};
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
    my $fh = $self->create_fh($self->filtered_mate1, ">>");
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
    my $fh = $self->create_fh($self->filtered_mate2, ">>") if($self->filtered_mate2);
    $self->{filtered_mate2_fh} = $fh;
  }
  return  $self->{filtered_mate2_fh};
}


sub close_files{
  my ($self) = @_;
  $self->{mate1_fh} = undef;
  $self->{mate2_fh} = undef;
  $self->{frag_fh} = undef;
  $self->{filtered_mate1_fh} = undef;
  $self->{filtered_mate2_fh} = undef;
  $self->{filtered_frag_fh} = undef;
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
    my $fh = $self->create_fh($self->filtered_frag, ">>") if($self->filtered_frag);
    $self->{filtered_frag_fh} = $fh;
  }
  return $self->{filtered_frag_fh};
}


sub create_filtered_filename{
  my ($self, $name) = @_;
  $name =~ s/\.fastq\.gz/\.filt\.fastq/;
  return $name;
}

sub create_fh{
  my ($self, $file, $mode, $use_pipe) = @_;
  my $fh = new FileHandle;
  my $open_cmd = $mode." ".$file;
  $open_cmd .= " | " if($use_pipe);
  $fh->open($open_cmd);
  if(!$fh){
    throw("Failed to open ".$open_cmd);
  }
  return $fh;
}


1;
