package ReseqTrack::Tools::RunAlignmentUtils;

use strict;
use warnings;

use ReseqTrack::Tools::FastQ qw (sample);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils;
use Data::Dumper;

use File::Basename;



use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
    get_base_read_counts
    subsample_fastq
    get_program_version  
    decide_file_skip     
);


sub get_base_read_counts {
  my $runalignment  = shift;
  my $collection = $runalignment->input;

  if  ( ! $collection->isa('ReseqTrack::Collection') ){
    print "Input is not a collection. Not getting base counts.\n";
    return;
  }
 
  my $files = $collection->others;

  my $bases = 0;
  my $reads = 1;
  my %base_counts;
  my %read_counts;
  my %read_lengths;
  my $max_read_length = 0;
  my $total_collection_bases = 0;
  my $total_collection_reads = 0;

  my $verbose =1 ;

  foreach my $f (@$files) {
    my $check = $f->statistics;
    foreach my $r (@$check) {
      if ( $r->{attribute_name} eq "base_count" ) {
	$bases = $r->{attribute_value};
	$total_collection_bases += $bases;
      }

      if ( $r->{attribute_name} eq "read_count" ) {
	$reads = $r->{attribute_value};
	$total_collection_reads += $reads;
      }

    }
    my $r_length= 0;
    $r_length = $bases / $reads if ($reads);

    if ( $r_length > $max_read_length) {
      $max_read_length = $r_length; 
    }

    print basename($f->name) ,"   $bases $reads $r_length\n";


    if ( $f->name eq $runalignment->fragment_file){
       $read_lengths{ $runalignment->fragment_file } = $bases / $reads if ($reads);
    }
    
    if (defined $runalignment->mate1_file){
      if ( $f->name eq $runalignment->mate1_file){
	$read_lengths{ $runalignment->mate1_file } = $bases / $reads if ($reads);
      }
    }

     if (defined $runalignment->mate2_file){
       if ( $f->name eq $runalignment->mate2_file){
	 $read_lengths{ $runalignment->mate2_file } = $bases / $reads if ($reads);
       }
     }

    $base_counts{ $f->name }  = $bases;
    $read_counts{ $f->name }  = $reads;
   # $read_lengths{ $f->name } = $bases / $reads if ($reads);
  }

  $runalignment->base_counts( \%base_counts );
  $runalignment->read_counts( \%read_counts );
  $runalignment->read_lengths( \%read_lengths );
  $runalignment->collection_base_count($total_collection_bases);

  print "Max read length = $max_read_length\n";
  $runalignment->max_read_length($max_read_length);


  return;
}






sub get_program_version{
  my ($runalignment)       = shift;
  my @aa;
  my $version;
 
#  print "Getting program version\n";

  $runalignment->program_version("UNK");


  if ( $runalignment->program ){
    my $cmd = $runalignment->program;
    @aa = `$cmd 2>&1`;
    chomp $aa[0] if @aa;
  }

  foreach my $i (@aa){
    if ( $i =~ /version/i){
      my @bb = split /Version\:/,$i ;
      $version = $bb[1];
      last;
    }
  }

  if ($version) {
    $version =~ s/\s+//;
    $version =~ s/\(.*\)//;
 #   print $version,"\n" if  $version;
    chomp $version;
    $runalignment->program_version($version); 
  }

#BFAST:   the blat-like fast accurate search tool
#Version: 0.6.4e git:Revision: undefined$

#Program: bwa (alignment via Burrows-Wheeler transformation)
#Version: 0.5.8c (r1536)

#shrimp2sam: LETTER SPACE (454,Illumina/Solexa,etc.) SPACE.
#SHRiMP 2.0.4 [ICC Intel(R) C++ g++ 4.1 mode]
-
return;
}

sub decide_file_skip{


  my ($runalignment)       = shift;
  $runalignment->skip_fragment("0");
  $runalignment->skip_mate_files("1");
  next if(!$runalignment->fragment_file || !(-e $runalignment->fragment_file));

  my $file_base_counts       = $runalignment->base_counts;
  my $collection_base_count  = $runalignment->collection_base_count;
 
  my $frag_bases = $$file_base_counts{$runalignment->fragment_file};  
  my $frag_percent_bases = ($frag_bases/$collection_base_count) *100;
  

  print "Fragment file bases \% collection bases =  $frag_percent_bases\n";
  if ($frag_percent_bases <= 50.0){
    print "Fragment file < 50\% of collection bases. Skipping alignment\n";
    $runalignment->skip_fragment("1");
    $runalignment->skip_mate_files("0");
  }
  
 return;

}

=pod
=head2  subsample_fastq

  Function  : subsample fastq files if required. If read lengths are different
  if paired end files, use the file with the longest read length to get indices
  for sampling. Then use these indices on the other mate file. If fastq are sampled
  then tmp files are created and these are assigned to fragment_file, mate1_file or
  mate2_file where appropriate.
  
  Exceptions: none
  Example   : $run_alignment->subsample_fastq;

=cut


sub subsample_fastq {
  my ($runalignment)       = shift;
  my $base_counts  = $runalignment->base_counts;
  my $read_counts  = $runalignment->read_counts;
  my $read_lengths = $runalignment->read_lengths;
  my $base;
  my $tmp_file;
  my $max = $runalignment->subsample_size;

  my $tmp_dir = $runalignment->working_dir;
  my $seq_list;

  my $total_reads   = 0;
  my $reads_used  = 0;



  if ( $runalignment->fragment_file && $runalignment->skip_fragment) {
    print "Skipping sampling check on fragment file\n";
  }

  if ( $runalignment->fragment_file && !$runalignment->skip_fragment) {

     if ( $$base_counts{ $runalignment->fragment_file } > $max ) {

      print "Sampling fragment file: ", $runalignment->fragment_file, "  ",
    $$base_counts{ $runalignment->fragment_file }, "\n";
 
      my $base = basename( $runalignment->fragment_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;

      $seq_list =
    sample( $runalignment->fragment_file, $tmp_file, $max );

      if ( scalar (keys %$seq_list) ) {
    print "Have indices\n";
    sample( $runalignment->fragment_file,
               $tmp_file, $max, $seq_list );
    
        $reads_used  += scalar ( keys  %$seq_list);
        $total_reads += $$read_counts{ $runalignment->fragment_file };
    print $reads_used,"\t",$total_reads,"\n";
    print "+reads used ", scalar ( keys  %$seq_list)," frag file \n";

        $runalignment->fragment_file($tmp_file);
        $runalignment->files_to_delete($tmp_file);

      }
      else{

       $reads_used       += $$read_counts{ $runalignment->fragment_file };
       $total_reads      += $$read_counts{ $runalignment->fragment_file };
       $runalignment->percent_reads_used ($reads_used/$total_reads*100 ) ;
        print  $runalignment->percent_reads_used,"\n";
     }

    }
    else {     
      print "No need to sample fragment file\n";
      $reads_used       += $$read_counts{ $runalignment->fragment_file };
      $total_reads      += $$read_counts{ $runalignment->fragment_file };
   }

     $runalignment->percent_reads_used ($reads_used/$total_reads*100 ) ;
     print  $runalignment->percent_reads_used,"\n";
  }

  return if  ($runalignment->skip_mate_files);

  #should have both
  if (  !($runalignment->mate1_file) &&  !($runalignment->mate2_file ) ) {     
    print "No mate files to subsample\n";
    return;
  }

  if (  ($runalignment->mate1_file) &&  !($runalignment->mate2_file ) ) {
      
    print "Have mate1 file but not mate2 file\n";
    print "Something wrong\n";
    throw ("Missing a mate file\n");
  }

  if (  ($runalignment->mate2_file) &&  !($runalignment->mate1_file ) ) {
      
    print "Have mate2 file but not mate1 file\n";
    print "Something wrong\n";
    throw ("Missing a mate file\n");
  }




  if ( ( $$base_counts{ $runalignment->mate1_file } < $max ) &&
       ( $$base_counts{ $runalignment->mate2_file } < $max ) ) {
    print "No need to subsample mate files\n";
    $runalignment->percent_reads_used (100 );
    print  "percent bases used ",$runalignment->percent_reads_used,"\n";
    return;
  }


  #Can get mate files where read lengths are different. Use file with longest
  #read length to get indices for subsampling. 



  if ( $$read_lengths{ $runalignment->mate1_file } !=
       $$read_lengths{ $runalignment->mate2_file } ) {

    print "OK. We appear to have unequal read lengths.\n";
    print "mate1 ", $$read_lengths{ $runalignment->mate1_file }, "\t", "mate2 ",
      $$read_lengths{ $runalignment->mate2_file }, "\n";

    if ( $$read_lengths{ $runalignment->mate1_file } >=
     $$read_lengths{ $runalignment->mate2_file } ) {
      print "Getting sampling indices from ", $runalignment->mate1_file, "\n";

#      print "Sampling ", $runalignment->mate1_file, "\n";
      $base     = basename( $runalignment->mate1_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;
      $seq_list =
    sample( $runalignment->mate1_file, $tmp_file, $max );
 
      if ( scalar (keys %$seq_list) ) {
    sample( $runalignment->mate1_file, $tmp_file, $max, $seq_list );
    $runalignment->files_to_delete($tmp_file);

        $reads_used    += scalar (keys %$seq_list);
    $total_reads   += $$read_counts{ $runalignment->mate1_file };
    print  "reads used  ", scalar (keys %$seq_list)," mate 1\n";

    $runalignment->mate1_file($tmp_file);
    $base     = basename( $runalignment->mate2_file );
    $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
    $tmp_file =~ s/\/\//\//;
    
    sample( $runalignment->mate2_file, $tmp_file, $max, $seq_list );

    $runalignment->files_to_delete($tmp_file);

    $reads_used     += scalar (keys %$seq_list);
    $total_reads    += $$read_counts{ $runalignment->mate2_file };
    print  "reads used ", scalar (keys %$seq_list)," mate 2\n";

    $runalignment->mate2_file($tmp_file);
    print $reads_used,"\t",$total_reads,"\n";
    $runalignment->percent_reads_used ( int ($reads_used/$total_reads*100) );
    print   "percent reads used ", $runalignment->percent_reads_used,"\n"; 
    return;
      }
      else{
    $runalignment->percent_reads_used ( 100) ;
    print  "No indices 100 percent reads used\n";
    return;
      }
    
    }


    else {
      print "Getting sampling indices from ", $runalignment->mate2_file, "\n";

 #     print "Sampling ", $runalignment->mate2_file, "\n";
      $base     = basename( $runalignment->mate2_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;
      $seq_list =
    sample( $runalignment->mate2_file, $tmp_file, $max );

      if ( scalar ( keys %$seq_list) ) {

    sample( $runalignment->mate2_file, $tmp_file, $max, $seq_list );
    

        $reads_used     += scalar keys(%$seq_list);
    $total_reads    += $$read_counts{ $runalignment->mate2_file };
    print  "reads $reads_used ", scalar keys(%$seq_list)," mate 2\n";

    $runalignment->mate2_file($tmp_file);
    $runalignment->files_to_delete($tmp_file);


    $base     = basename( $runalignment->mate1_file );
    $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
    $tmp_file =~ s/\/\//\//;
    sample( $runalignment->mate1_file, $tmp_file, $max, $seq_list );
    print  "reads $reads_used ", scalar keys(%$seq_list)," mate 1\n";

    $reads_used = scalar ( keys %$seq_list);
    $total_reads    += $$read_counts{ $runalignment->mate1_file };

    $runalignment->mate1_file($tmp_file);
    $runalignment->files_to_delete($tmp_file);

        print $reads_used,"\t",$total_reads,"\n";
    $runalignment->percent_reads_used ( int ($reads_used/$total_reads*100) );
    print $reads_used,"\t",$total_reads,"\n";
    print   "percent reads used ", $runalignment->percent_reads_used,"\n";
    return;
      }
      else{
    $runalignment->percent_reads_used ( 100) ;
    print  "No indices 100 percent reads used\n";
      }
    }
    return;
  }
    
 


 # RES. I think below is overkill.

  
  if ( $$base_counts{ $runalignment->mate1_file } >=  $$base_counts{ $runalignment->mate2_file } ) {
    print "Sampling: ", $runalignment->mate1_file, "\n";
    $base = basename( $runalignment->mate1_file );
    $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
    $tmp_file =~ s/\/\//\//;  
    $seq_list =
      sample( $runalignment->mate1_file, $tmp_file, $max ) ;

    if (scalar (keys %$seq_list) ) {
      print "\nUsing indices from for sampling :\n" . $runalignment->mate1_file. "\n";
      sample( $runalignment->mate1_file, $tmp_file, $max, $seq_list );
    
 #     print "\nself mate1 ",$runalignment->mate1_file,"\n";
      $reads_used     += scalar  (keys %$seq_list);
      $total_reads    +=  $$read_counts{ $runalignment->mate1_file };
      print  "reads $reads_used ",  scalar  (keys %$seq_list)," mate 1\n";
      
      $runalignment->mate1_file($tmp_file);
      $runalignment->files_to_delete($tmp_file);


      $base = basename( $runalignment->mate2_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;
      sample( $runalignment->mate2_file, $tmp_file, $max, $seq_list );    
      $reads_used = scalar  (keys %$seq_list);
     
      $reads_used     += scalar  (keys %$seq_list);
      $total_reads    +=  $$read_counts{ $runalignment->mate2_file };
      print  "reads $reads_used ",  scalar  (keys %$seq_list)," mate 2\n";
       
      $runalignment->mate2_file($tmp_file);
      $runalignment->files_to_delete($tmp_file); 

      print $reads_used,"\t",$total_reads,"\n";
      $runalignment->percent_reads_used ( int ($reads_used/$total_reads*100) );
      print  "percent reads used ", $runalignment->percent_reads_used,"\n";
    }
    else{
      $runalignment->percent_reads_used ( 100) ;
      print  "No indices 100 percent reads used\n";
    }
      return;
  }



   
  if ( $$base_counts{ $runalignment->mate2_file } >  $$base_counts{ $runalignment->mate2_file } ) {
    print "Sampling: ", $runalignment->mate2_file, "\n";
    $base = basename( $runalignment->mate2_file );
    $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
    $tmp_file =~ s/\/\//\//;


    $seq_list =
      sample( $runalignment->mate2_file, $tmp_file, $max ) ;

    if (scalar (keys %$seq_list) ) {
      print "Using indices from for sampling:\n" . $runalignment->mate2_file. "\n";
    
      sample( $runalignment->mate2_file, $tmp_file, $max, $seq_list );

      $reads_used     += scalar  (keys %$seq_list);
      $total_reads    +=  $$read_counts{ $runalignment->mate2_file };
      print  "reads $reads_used ",  scalar  (keys %$seq_list)," mate 2\n";

      $runalignment->mate2_file($tmp_file);
      $runalignment->files_to_delete($tmp_file);

      $base = basename( $runalignment->mate1_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;
 #     print "$tmp_file\n";
      sample( $runalignment->mate1_file, $tmp_file, $max, $seq_list );
    
      $reads_used     += scalar  (keys %$seq_list);
      $total_reads    +=  $$read_counts{ $runalignment->mate1_file };
      print  "reads $reads_used ",  scalar  (keys %$seq_list)," mate 2\n";
      $runalignment->mate1_file($tmp_file);
      $runalignment->files_to_delete($tmp_file);
      print $reads_used,"\t",$total_reads,"\n";
      $runalignment->percent_reads_used (int ($reads_used/$total_reads*100) );
      print   "percent reads used ",$runalignment->percent_reads_used," ++++\n";
    } 
    else{
      $runalignment->percent_reads_used ( 100) ;
      print  "No indices 100 percent reads used\n";
    }
    return;

 }

  return;
}



1;
