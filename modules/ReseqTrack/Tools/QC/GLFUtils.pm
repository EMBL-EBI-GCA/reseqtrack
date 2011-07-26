package ReseqTrack::Tools::QC::GLFUtils;

use strict;
use warnings;
use Exporter;

use Data::Dumper;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::GenotypeResultsAdaptor;


use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::SequenceIndexUtils qw (assign_files);
use ReseqTrack::Tools::FileUtils qw (get_count_stats);
use ReseqTrack::Tools::FastqUtils  qw ( sample ); 
use File::Basename;
use Data::Dumper;
use vars qw (@ISA  @EXPORT_OK);

@ISA       = qw(Exporter);
@EXPORT_OK = qw(
  check_previous_result
  skip_lane
  decide_file_skip
  have_snps
  get_max_read_length 
  sample_array_list_of_files
  skip_lane
  get_bam_mapped_stats
);

sub check_previous_result{
  my ($gr, $name,$update) = @_;

  
  my $prev_result = $gr->fetch_by_name ($name);
  if ($prev_result) {  
    throw ( "Already have results for $name. 'update' option not used")
      unless $update;
  } else {
    print "No current results in DB\n";
  }
  return $prev_result;
} 


sub skip_lane{

  my ( $db,$other_id,$name,$claimed_sample,$prev_result,$reason,$input) = @_;
  
  print "SKIPPING LANE\n";
  # print " $run_alignment,$db,$other_id,$name,$claimed_sample,$prev_result,$reason\n";

  
  my $GTR ="ReseqTrack::GenotypeResults";
  my $genotype_results = $GTR->new(
				   -other_id      => $other_id,
				   -table_name    => "collection",
				   -name          => $name,
				   -claimed       => $claimed_sample,
				   -top_hit       => "N/A",
				   -second_hit    => "N/A",
				   -ratio21       => "N/A",
				   -ratio_claimed => "N/A",
				   -reference     => "N/A",
				   -snps_bin      => "N/A",
				   -aligner       => "N/A",
				   -validation_method =>  "N/A",
				   -percent_mapped => 0,
				   -percent_reads_used => 0,
				   -max_bases     => 0,
				   -verdict       => $reason, 
				   -cfg_file      => "N/A",
				   -skip_others   => '1',
				  );
	
	

  if (! $db) {
    $db = ReseqTrack::DBSQL::DBAdaptor->new(
					    -host   => $$input{dbhost},
					    -user   => $$input{dbuser},
					    -port   => $$input{dbport},
					    -dbname => $$input{dbname},
					    -pass   => $$input{dbpass},
					   );	
  }
	
  my  $gra = $db->get_GenotypeResultsAdaptor;



  if ($prev_result) {
    print "Updating\n";
    $gra->update($genotype_results) ;
  } else {
    print "Storing\n";
    $gra->store($genotype_results) ;
  }
  
  print "\nSkipping $name:$reason run\n";
  
  return 1;
}



sub  decide_file_skip {

  my ($collection, $subsample_limit)  = @_;

  if ( ! $collection->isa('ReseqTrack::Collection') ) {
    print "Input is not a collection. Not getting base counts.\n";
    return;
  }

  my $files = $collection->others;

  my $bases = 0;
  my $reads = 1;
  my %base_counts;
  my $total_collection_bases = 0;
  my $frag_per = 0.0;
  my $mate_per = 0.0;
  my $verbose =1 ;
  my @tmp;
  my $align_these;
  my  ($mate1, $mate2, $frag);
  $mate1 = "";
  $mate2 = "";
  $frag = "";

  my $requires_subsampling = 0;


  foreach my $f (@$files) {
    push (@tmp, $f->name);
  }



  ( $mate1, $mate2, $frag ) = assign_files( \@tmp );

  foreach my $f (@$files) {
    my  ($read_count, $base_count) = get_count_stats ($f);
    $base_counts{ $f->name }  =  $base_count;
    $total_collection_bases  +=  $base_count;
  }
  
  
  if ( defined  $base_counts{$frag}) {
    print "Frag file basecount \% ";
    $frag_per =  $base_counts{$frag} / $total_collection_bases * 100;
    print "$frag_per \n";
  }

  if ($mate1 && $mate2 ) {
    print "Mate file basecount \% ";
    $mate_per =  ($base_counts{$mate1}+ $base_counts{$mate2}) / $total_collection_bases * 100;
    print "$mate_per\n";
  }
  
  if ( $frag_per < 50) {
    push (@$align_these,  $mate1);
    push (@$align_these, $mate2);
  } else {
    push (@$align_these, $frag);
  }

  print "Aligning\n";
  foreach (@$align_these) {
    print "-- $_\n";

   $requires_subsampling++ if ($base_counts{$_} >  $subsample_limit) ;
  }

  print "No need to sample files\n" if ( ! $requires_subsampling);
  print "==============\n";
  return $align_these,  $requires_subsampling;  
}





sub have_snps{

  my ($sample, $snps_list) = @_;

  if ($sample =~ /unidentified/) {
    print "sample unidentified";
    sleep (2);
    throw ("sample unidentified");
  }

  throw ("Could not find snp list") unless (-e $snps_list);

  my @aa = `grep $sample $snps_list`;
  chomp $aa[0] if @aa;

  if (@aa) {
    print "Have snps for $sample\n";
    return 1;
  }

  print"\nDo not have snps for $sample\n";
  

  return 0;
}





sub  get_max_read_length {

  my ($collection, $aligning, $verbose )   = @_;
  my $max_read_length = 0;
  my @using; 

  if ( ! $collection->isa('ReseqTrack::Collection') ) {
    print "Input is not a collection. Not getting base counts.\n";
    return;
  }

  my $files = $collection->others;

  foreach my $f (@$files) {
    print  $f->name,"\n";
    my $name = $f->name;
    @using =  grep (/$name/, @$aligning);
    next if (! @using);

    my  ($read_count, $base_count) = get_count_stats ($f);

    print $f->filename, " reads $read_count bases $base_count\n" if $verbose;;

    my $this_read_length = $base_count/$read_count;
    if ( $this_read_length > $max_read_length) {
      $max_read_length =  $this_read_length;
    }
  }
  
  print "max_read_length $max_read_length\n";
  
  die "Could not calc max readlength\n" if ($max_read_length == 0);

  return $max_read_length;
}




sub sample_array_list_of_files{

  my ( $files, $max, $tmp_dir) = @_;
  my @files_to_delete;
  my @files_to_process;
  my ($tmp_file, $base);
  my $seq_list;

  my ( $mate1, $mate2, $fragment_file ) = assign_files ( $files);

  if ( $fragment_file ) {
    print "Sampling fragment file: ", $fragment_file,  "\n";
 
    $base = basename( $fragment_file );
    $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
    $tmp_file =~ s/\/\//\//;

    $seq_list = sample( $fragment_file, $tmp_file, $max );

    if ( scalar (keys %$seq_list) ) {
      print "Have indices\n";
      sample( $fragment_file, $tmp_file, $max, $seq_list );
    
      push (@files_to_process ,$tmp_file);
      push (@files_to_delete  ,$tmp_file);
    }

  }

  #should have both
  if (  !($mate1) &&  !($mate2 ) ) {     
    print "No mate files to subsample\n";
    return ( \@files_to_process, \@files_to_delete)  ;
  }

  if (  ($mate1) &&  !($mate2 ) ) {
      
    print "Have mate1 file but not mate2 file\n";
    print "Something wrong\n";
    throw ("Missing a mate file\n");
  }

  if (  ($mate2) &&  !($mate1 ) ) {
      
    print "Have mate2 file but not mate1 file\n";
    print "Something wrong\n";
    throw ("Missing a mate file\n");
  }


 
  $base     = basename( $mate1);
  $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
  $tmp_file =~ s/\/\//\//;
  $seq_list =	sample( $mate1, $tmp_file, $max );
 
  if ( scalar (keys %$seq_list) ) {
    sample( $mate1, $tmp_file, $max, $seq_list );
    push (@files_to_process ,$tmp_file);
    push (@files_to_delete  ,$tmp_file);
     
    $base     = basename( $mate2 );
    $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
    $tmp_file =~ s/\/\//\//;   
    sample( $mate2, $tmp_file, $max, $seq_list );     
    push (@files_to_process ,$tmp_file);
    push (@files_to_delete  ,$tmp_file);
  }    


  if (!@files_to_process) {
    return $files;
  }

  return ( \@files_to_process, \@files_to_delete)  ;


}


sub get_bam_mapped_stats {

  my ($samtools, $bam)       = @_;
   
  my $flagstat =  $samtools . " flagstat ";

  my $mapped = 0;

  my $cmd = $flagstat . " " .$bam;
  my @stats = `$cmd`;

 

  foreach (@stats) {
    print;
    if ( $_ =~ /mapped \(/) {
      my @aa = split /\s+/;
      $mapped = $aa[-1];
      $mapped =~ s/\(|\)|\%//g;
    }
  }
        
  return $mapped;
}

1;
