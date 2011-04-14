
=pod

=head1 NAME

ReseqTrack::Tools::RunAlignment

=head1 SYNOPSIS

This is a base class for RunAlignment objects and provides some standard
accessor methods and throws exceptions when vital methods aren't implemented in
the child classes. The Child classes should wrap specific alignment algorithms

=head1 Example


=cut

package ReseqTrack::Tools::RunAlignment;

use strict;
use warnings;

use ReseqTrack::Tools::FastQ qw (sample);
use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::Argument qw(rearrange);
use ReseqTrack::Tools::SequenceIndexUtils;
use Data::Dumper;

use File::Basename;


=head2 new

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, path to reference genome file/directory/stem as is required
  for aligner
  Arg [3]   : string/arrayref/ReseqTrack::File/ReseqTrack::Collection. The
  input argument can be several things, all must provide a filename or set of 
  filenames for the aligner to run on
  Arg [4]   : string, path to alignment program
  Arg [5]   : string, commandline options for program
  Arg [6]   : path to samtools installation so output bam can be
  indexed
  Function  : This should create the ReseqTrack::Tools::RunAlignment object 
  Returntype: ReseqTrack::Tools::RunAlignment
  Exceptions: 
  Example   : 

=cut

sub new {
  my ( $class, @args ) = @_;
  my $self = {};
  bless $self, $class;

  #reference genome file
  #input sequence, we will allow single file name, list of file names
  # or a collection object with associated sequence
  #program
  #commandline options
  my (
      $reference, $input,     $program,
      $options,   $samtools,  $working_dir,
      $name,      $file_info, $subsample_size,
      $run_meta_info)
    =

      rearrange(
		[
		 qw(
              REFERENCE
              INPUT
              PROGRAM
              OPTIONS
              SAMTOOLS
              WORKING_DIR
              NAME
              FILE_INFO
              SUBSAMPLE_SIZE
              RUN_META_INFO
)
		],
		@args
	       );
  #some defaults
#  $self->subsample_size("1200000000");
  $self->working_dir("/tmp/") unless ($working_dir);
  $self->fragment_file ("");
  $self->mate1_file("");
  $self->mate2_file("");

  #####


  $self->reference($reference);
  $self->input($input);
  $self->program($program);
  $self->options($options);
  $self->samtools($samtools);
  $self->working_dir($working_dir);


  $self->subsample_size($subsample_size);
  $self->file_info($file_info);
  $self->name($name);
  $self->run_meta_info($run_meta_info);


  unless ( $self->name ) {
    my $string = $self->mate1_file;
    $string = $self->fragment_file unless ($string);
    $string =~ /^(\S+)\./;
    $self->name($1);
    throw(
	  "ReseqTrack::Tools::RunAlignment, Not sure what to do failed to define "
	  . "a name from "
	  . $string
	  . " for this run" )
      unless ( $self->name );
  }

  #
  throw("Have no working directory")
    unless ( $self->working_dir && -d $self->working_dir );


  $self->get_base_read_counts();

  $self->get_program_version();


  
  return $self;
}
=pod
=head2


=cut

sub get_program_version{
  my ($self)       = shift;
  my @aa;
  my $version;
 
#  print "Getting program version\n";

  $self->program_version("UNK");


  if ( $self->program ){
    my $cmd = $self->program;
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
    $self->program_version($version); 
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


  my ($self)       = shift;

  my $file_base_counts       = $self->base_counts;
  my $collection_base_count  = $self->collection_base_count;
 
  my $frag_bases = $$file_base_counts{$self->fragment_file};  
  my $frag_percent_bases = ($frag_bases/$collection_base_count) *100;
  

  print "Fragment file bases \% collection bases =  $frag_percent_bases\n";
  if ($frag_percent_bases <= 50.0){
    print "Fragment file < 50\% of collection bases. Skipping alignment\n";
    $self->skip_fragment("1");
    $self->skip_mate_files("0");
  }
  else{
    print "Mate files < 50\% of collection bases. Skipping alignment\n";
    $self->skip_fragment("0");
    $self->skip_mate_files("1");

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
  my ($self)       = shift;
  my $base_counts  = $self->base_counts;
  my $read_counts  = $self->read_counts;
  my $read_lengths = $self->read_lengths;
  my $base;
  my $tmp_file;
  my $max = $self->subsample_size;

  my $tmp_dir = $self->working_dir;
  my $seq_list;

  my $total_reads   = 0;
  my $reads_used  = 0;



  if ( $self->fragment_file && $self->skip_fragment) {
    print "Skipping sampling check on fragment file\n";
  }

  if ( $self->fragment_file && !$self->skip_fragment) {

     if ( $$base_counts{ $self->fragment_file } > $max ) {

      print "Sampling fragment file: ", $self->fragment_file, "  ",
	$$base_counts{ $self->fragment_file }, "\n";
 
      my $base = basename( $self->fragment_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;

      $seq_list =
	sample( $self->fragment_file, $tmp_file, $max );

      if ( scalar (keys %$seq_list) ) {
	print "Have indices\n";
	sample( $self->fragment_file,
		       $tmp_file, $max, $seq_list );
	
        $reads_used  += scalar ( keys  %$seq_list);
        $total_reads += $$read_counts{ $self->fragment_file };
	print $reads_used,"\t",$total_reads,"\n";
	print "+reads used ", scalar ( keys  %$seq_list)," frag file \n";

        $self->fragment_file($tmp_file);
        $self->files_to_delete($tmp_file);

      }
      else{

       $reads_used       += $$read_counts{ $self->fragment_file };
       $total_reads      += $$read_counts{ $self->fragment_file };
       $self->percent_reads_used ($reads_used/$total_reads*100 ) ;
        print  $self->percent_reads_used,"\n";
     }

    }
    else {     
      print "No need to sample fragment file\n";
      $reads_used       += $$read_counts{ $self->fragment_file };
      $total_reads      += $$read_counts{ $self->fragment_file };
   }

     $self->percent_reads_used ($reads_used/$total_reads*100 ) ;
     print  $self->percent_reads_used,"\n";
  }

  return if  ($self->skip_mate_files);

  #should have both
  if (  !($self->mate1_file) &&  !($self->mate2_file ) ) {     
    print "No mate files to subsample\n";
    return;
  }

  if (  ($self->mate1_file) &&  !($self->mate2_file ) ) {
      
    print "Have mate1 file but not mate2 file\n";
    print "Something wrong\n";
    throw ("Missing a mate file\n");
  }

  if (  ($self->mate2_file) &&  !($self->mate1_file ) ) {
      
    print "Have mate2 file but not mate1 file\n";
    print "Something wrong\n";
    throw ("Missing a mate file\n");
  }




  if ( ( $$base_counts{ $self->mate1_file } < $max ) &&
       ( $$base_counts{ $self->mate2_file } < $max ) ) {
    print "No need to subsample mate files\n";
    $self->percent_reads_used (100 );
    print  "percent bases used ",$self->percent_reads_used,"\n";
    return;
  }


  #Can get mate files where read lengths are different. Use file with longest
  #read length to get indices for subsampling. 



  if ( $$read_lengths{ $self->mate1_file } !=
       $$read_lengths{ $self->mate2_file } ) {

    print "OK. We appear to have unequal read lengths.\n";
    print "mate1 ", $$read_lengths{ $self->mate1_file }, "\t", "mate2 ",
      $$read_lengths{ $self->mate2_file }, "\n";

    if ( $$read_lengths{ $self->mate1_file } >=
	 $$read_lengths{ $self->mate2_file } ) {
      print "Getting sampling indices from ", $self->mate1_file, "\n";

#      print "Sampling ", $self->mate1_file, "\n";
      $base     = basename( $self->mate1_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;
      $seq_list =
	sample( $self->mate1_file, $tmp_file, $max );
 
      if ( scalar (keys %$seq_list) ) {
	sample( $self->mate1_file, $tmp_file, $max, $seq_list );
	$self->files_to_delete($tmp_file);

        $reads_used    += scalar (keys %$seq_list);
	$total_reads   += $$read_counts{ $self->mate1_file };
	print  "reads used  ", scalar (keys %$seq_list)," mate 1\n";

	$self->mate1_file($tmp_file);
	$base     = basename( $self->mate2_file );
	$tmp_file = $tmp_dir . "/" . "$$\_" . $base;
	$tmp_file =~ s/\/\//\//;
	
	sample( $self->mate2_file, $tmp_file, $max, $seq_list );

	$self->files_to_delete($tmp_file);

	$reads_used     += scalar (keys %$seq_list);
	$total_reads    += $$read_counts{ $self->mate2_file };
	print  "reads used ", scalar (keys %$seq_list)," mate 2\n";

	$self->mate2_file($tmp_file);
	print $reads_used,"\t",$total_reads,"\n";
	$self->percent_reads_used ( int ($reads_used/$total_reads*100) );
	print   "percent reads used ", $self->percent_reads_used,"\n"; 
	return;
      }
      else{
	$self->percent_reads_used ( 100) ;
	print  "No indices 100 percent reads used\n";
	return;
      }
    
    }


    else {
      print "Getting sampling indices from ", $self->mate2_file, "\n";

 #     print "Sampling ", $self->mate2_file, "\n";
      $base     = basename( $self->mate2_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;
      $seq_list =
	sample( $self->mate2_file, $tmp_file, $max );

      if ( scalar ( keys %$seq_list) ) {

	sample( $self->mate2_file, $tmp_file, $max, $seq_list );
	

        $reads_used     += scalar keys(%$seq_list);
	$total_reads    += $$read_counts{ $self->mate2_file };
	print  "reads $reads_used ", scalar keys(%$seq_list)," mate 2\n";

	$self->mate2_file($tmp_file);
	$self->files_to_delete($tmp_file);


	$base     = basename( $self->mate1_file );
	$tmp_file = $tmp_dir . "/" . "$$\_" . $base;
	$tmp_file =~ s/\/\//\//;
	sample( $self->mate1_file, $tmp_file, $max, $seq_list );
	print  "reads $reads_used ", scalar keys(%$seq_list)," mate 1\n";

	$reads_used = scalar ( keys %$seq_list);
	$total_reads    += $$read_counts{ $self->mate1_file };

	$self->mate1_file($tmp_file);
	$self->files_to_delete($tmp_file);

      	print $reads_used,"\t",$total_reads,"\n";
	$self->percent_reads_used ( int ($reads_used/$total_reads*100) );
	print $reads_used,"\t",$total_reads,"\n";
	print   "percent reads used ", $self->percent_reads_used,"\n";
	return;
      }
      else{
	$self->percent_reads_used ( 100) ;
	print  "No indices 100 percent reads used\n";
      }
    }
	return;
  }
    
 


 # RES. I think below is overkill.

  
  if ( $$base_counts{ $self->mate1_file } >=  $$base_counts{ $self->mate2_file } ) {
    print "Sampling: ", $self->mate1_file, "\n";
    $base = basename( $self->mate1_file );
    $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
    $tmp_file =~ s/\/\//\//;  
    $seq_list =
      sample( $self->mate1_file, $tmp_file, $max ) ;

    if (scalar (keys %$seq_list) ) {
      print "\nUsing indices from for sampling :\n" . $self->mate1_file. "\n";
      sample( $self->mate1_file, $tmp_file, $max, $seq_list );
    
 #     print "\nself mate1 ",$self->mate1_file,"\n";
      $reads_used     += scalar  (keys %$seq_list);
      $total_reads    +=  $$read_counts{ $self->mate1_file };
      print  "reads $reads_used ",  scalar  (keys %$seq_list)," mate 1\n";
      
      $self->mate1_file($tmp_file);
      $self->files_to_delete($tmp_file);


      $base = basename( $self->mate2_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;
      sample( $self->mate2_file, $tmp_file, $max, $seq_list );    
      $reads_used = scalar  (keys %$seq_list);
     
      $reads_used     += scalar  (keys %$seq_list);
      $total_reads    +=  $$read_counts{ $self->mate2_file };
      print  "reads $reads_used ",  scalar  (keys %$seq_list)," mate 2\n";
       
      $self->mate2_file($tmp_file);
      $self->files_to_delete($tmp_file); 

      print $reads_used,"\t",$total_reads,"\n";
      $self->percent_reads_used ( int ($reads_used/$total_reads*100) );
      print  "percent reads used ", $self->percent_reads_used,"\n";
    }
    else{
      $self->percent_reads_used ( 100) ;
      print  "No indices 100 percent reads used\n";
    }
      return;
  }



   
  if ( $$base_counts{ $self->mate2_file } >  $$base_counts{ $self->mate2_file } ) {
    print "Sampling: ", $self->mate2_file, "\n";
    $base = basename( $self->mate2_file );
    $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
    $tmp_file =~ s/\/\//\//;


    $seq_list =
      sample( $self->mate2_file, $tmp_file, $max ) ;

    if (scalar (keys %$seq_list) ) {
      print "Using indices from for sampling:\n" . $self->mate2_file. "\n";
    
      sample( $self->mate2_file, $tmp_file, $max, $seq_list );

      $reads_used     += scalar  (keys %$seq_list);
      $total_reads    +=  $$read_counts{ $self->mate2_file };
      print  "reads $reads_used ",  scalar  (keys %$seq_list)," mate 2\n";

      $self->mate2_file($tmp_file);
      $self->files_to_delete($tmp_file);

      $base = basename( $self->mate1_file );
      $tmp_file = $tmp_dir . "/" . "$$\_" . $base;
      $tmp_file =~ s/\/\//\//;
 #     print "$tmp_file\n";
      sample( $self->mate1_file, $tmp_file, $max, $seq_list );
    
      $reads_used     += scalar  (keys %$seq_list);
      $total_reads    +=  $$read_counts{ $self->mate1_file };
      print  "reads $reads_used ",  scalar  (keys %$seq_list)," mate 2\n";
      $self->mate1_file($tmp_file);
      $self->files_to_delete($tmp_file);
      print $reads_used,"\t",$total_reads,"\n";
      $self->percent_reads_used (int ($reads_used/$total_reads*100) );
      print   "percent reads used ",$self->percent_reads_used," ++++\n";
    } 
    else{
      $self->percent_reads_used ( 100) ;
      print  "No indices 100 percent reads used\n";
    }
    return;

 }

  return;
}

sub get_base_read_counts {
  my $self  = shift;
  my $files = $self->file_info;
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


    if ( $f->name eq $self->fragment_file){
       $read_lengths{ $self->fragment_file } = $bases / $reads if ($reads);
    }
    
    if (defined $self->mate1_file){
      if ( $f->name eq $self->mate1_file){
	$read_lengths{ $self->mate1_file } = $bases / $reads if ($reads);
      }
    }

     if (defined $self->mate2_file){
       if ( $f->name eq $self->mate2_file){
	 $read_lengths{ $self->mate2_file } = $bases / $reads if ($reads);
       }
     }

    $base_counts{ $f->name }  = $bases;
    $read_counts{ $f->name }  = $reads;
   # $read_lengths{ $f->name } = $bases / $reads if ($reads);
  }

  $self->base_counts( \%base_counts );
  $self->read_counts( \%read_counts );
  $self->read_lengths( \%read_lengths );
  $self->collection_base_count($total_collection_bases);

  print "Max read length = $max_read_length\n";
  $self->max_read_length($max_read_length);


  return;
}

sub base_counts {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'base_counts'} = $arg;
  }
  return $self->{'base_counts'};
}

sub read_counts {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'read_counts'} = $arg;
  }
  return $self->{'read_counts'};
}

sub read_lengths {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'read_lengths'} = $arg;
  }
  return $self->{'read_lengths'};
}

sub file_info {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'file_info'} = $arg;
  } else {

  }
  return $self->{'file_info'};
}

=head2 accessor methods

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, generally
  Function  : These are accessor methods for variables needed to run
  and alignment like reference genomes, program and options
 Returntype: string, generally
 Exceptions: n/a
  Example   : my $reference = $run_alignment->reference;

=cut

sub reference {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'reference'} = $arg;
  }
  return $self->{'reference'};
}

sub program {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'program'} = $arg;
  }
  return $self->{'program'};
}

sub options {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'options'} = $arg;
  }
  return $self->{'options'};
}

sub samtools {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'samtools'} = $arg;
  }
  return $self->{'samtools'};
}

sub fragment_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'fragment_file'} = $arg;
  }
  return $self->{'fragment_file'};
}

sub mate1_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'mate1_file'} = $arg;
  }
  return $self->{'mate1_file'};
}

sub mate2_file {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'mate2_file'} = $arg;
  }
  return $self->{'mate2_file'};
}

sub name {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'name'} = $arg;
  }
  return $self->{'name'};
}
sub max_read_length{
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'max_read_length'} = $arg;
  }
  return $self->{'max_read_length'};

}
sub subsample_size {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{'subsample_size'} = $arg;
  }
  return $self->{'subsample_size'};
}

sub files_to_delete {
  my ( $self, $file ) = @_;
  if ($file) {
    if ( ref($file) eq 'ARRAY' ) {
      foreach my $path (@$file) {
	$self->{'files_to_delete'}->{$path} = 1;
      }
    } else {
      $self->{'files_to_delete'}->{$file} = 1;
    }
  }
  my @keys = keys( %{ $self->{'files_to_delete'} } );
  return \@keys;
}

=head2 working_dir

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, path of working directory
  Function  : accessor method for working directory, also creates directory
  if it doesn't exist
  Returntype: string
  Exceptions: n/a
  Example   : my $work_dir = $self->working_dir;

=cut

sub working_dir {
  my ( $self, $arg ) = @_;

  if ($arg) {
    $self->{'working_dir'} = $arg;
  }

  if ( $self->{'working_dir'} ) {
    unless ( -d $self->{'working_dir'} ) {
      mkdir( $self->{'working_dir'}, 775 );
    }
    unless ( -d $self->{'working_dir'} ) {
      throw(  "ReseqTrack::Tools::RunAlignment cannot run when "
	      . $self->{'working_dir'}
	      . " does not exist " );
    }
  }
  return $self->{'working_dir'};
}

=head2 change_dir

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, path of working directory
  Function  : changes to given working directory
  Returntype: string
  Exceptions: throws if it can't change to given directory
  Example   : $self->check_dir();

=cut

sub change_dir {
  my ( $self, $dir ) = @_;
  $dir = $self->working_dir unless ($dir);
  chdir($dir)
    or throw( "Failed to change to " 
	      . $dir
	      . " ReseqTrack::Tools::RunAlignment check_dir" );
  return $dir;
}

=head2 output_files

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string or arrayref of strings
  Function  : This is to store output filepaths
  Returntype: arrayref of strings
  Exceptions: n/a
  Example   : $self->output('path/to/file');

=cut

sub output_files {
  my ( $self, $arg ) = @_;
  $self->{'output'} = [] unless ( $self->{'output'} );
  if ($arg) {
    if ( ref($arg) eq 'ARRAY' ) {
      push( @{ $self->{'output'} }, @$arg );
    } else {
      push( @{ $self->{'output'} }, $arg );
    }
  }
  return $self->{'output'};
}

=head2 input

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string/arrayref/ReseqTrack::File/ReseqTrack::Collection. The
  input argument can be several things, all must provide a filename or set of 
  filenames for the aligner to run on
  Function  : setup the input file paths
  Returntype: string/arrayref/ReseqTrack::File/ReseqTrack::Collection, it returns
    what is passed in
  Exceptions: 
  Example   : 

=cut

sub input {
  my ( $self, $input ) = @_;

  #need to work out if we have filepaths, file objects or a collection object
  if ($input) {
    if ( -e $input ) {
      $self->fragment_file($input);
    } elsif ( ref($input) eq 'ARRAY' ) {
      my ( $mate1, $mate2, $frag );
      if ( -e $input->[0] ) {
	( $mate1, $mate2, $frag ) = $self->assign_fastq_files($input);
      } elsif ( $input->[0]->isa("ReseqTrack::File") ) {
	my @names;
	foreach my $file (@$input) {
	  push( @names, $file->name );
	}
	( $mate1, $mate2, $frag ) =
	  $self->assign_fastq_files( \@names );
      } else {
	print STDERR "Not sure how to deal with the contents of "
	  . $input . "\n";
	foreach my $element (@$input) {
	  print STDERR $element . "\n";
	}
	throw(
	      "ReseqTrack::Tools::RunAlignment::input Failed to process "
	      . $input );
      }
      $self->fragment_file($frag);
      $self->mate1_file($mate1);
      $self->mate2_file($mate2);
    } elsif ( $input->isa("ReseqTrack::File") ) {
      $self->fragment_file( $input->name );
    } elsif ( $input->isa('ReseqTrack::Collection') ) {
      throw(
	    "ReseqTrack::Tools::RunAlignment::input Can only handle file "
	    . "collections not "
	    . $input->table_name
	    . " collections" )
	unless ( $input->table_name eq 'file' );
      my $others = $input->others;
      my @names;
      foreach my $other (@$others) {
	push( @names, $other->name );
      }
      my ( $mate1, $mate2, $frag ) = $self->assign_fastq_files( \@names );
      $self->fragment_file($frag);
      $self->mate1_file($mate1);
      $self->mate2_file($mate2);
    } else {
      throw(
	    "ReseqTrack::Tools::RunAlignment::input not sure how to handle input "
	    . $input
	    . " what type is it?" );
    }
    $self->{'input'} = $input;
  }
  return $self->{'input'};
}

=head2 run

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Function  : each child object should implement a run method
  Returntype: n/a
  Exceptions: throws as this method should be implement in the child class
  Example   : 

=cut

sub run {
  my ($self) = @_;
  throw(  $self
          . " must implement a run method as ReseqTrack::Tools::RunAlignment "
          . "does not provide one" );
}

=head2 create_bam_from_sam

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : string, path to sam file
  Arg [3]   : binary flag if set sam files marked for deletion, this is on by
  default but it can be turned off if delete_sam is set to 0
  Function  : create a bam file from a sam file
  Returntype: string, path to bam file
  Exceptions: n/a
  Example   : 

=cut

sub create_bam_from_sam {
  my ( $self, $sam, $delete_sam ) = @_;
  unless ( defined($delete_sam) ) {
    $delete_sam = 1;
  }
  my $bam = $sam;
  $bam =~ s/sam/bam/;
  my $cmd =
    $self->samtools . " import " . $self->reference . " " . $sam . " " . $bam;
  print $cmd. "\n";
  eval {
    my $exit = system($cmd);
    if ( $exit && $exit >= 1 ) {
      throw( "Failed to run " . $cmd );
    }
  };
  if ($@) {
    throw("Failed to run samtools import $@");
  }
  if ($delete_sam) {
    $self->files_to_delete($sam);
  }
  return $bam;
}

=head2 delete_files

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : arrayref of string of paths to delete
  Function  : remove given files
  Returntype: n/a
  Exceptions: throws if a file doesn't exist
  Example   : 

=cut

sub delete_files {
  my ( $self, $files ) = @_;
  $files = $self->files_to_delete unless ($files);
  foreach my $file (@$files) {
    print "Deleting " . $file . "\n";
    unlink $file;
  }
  return;
}

=head2 assign_fastq_files

  Arg [1]   : ReseqTrack::Tools::RunAlignment
  Arg [2]   : arrayref of filepaths
  Function  : assign file names into mate1, mate2 and frag
  Returntype: array of filepaths
  Exceptions: if the path doesn't match any of the regexs
  Example   : my ($mate1, $mate2, $frag) = $self->assign_files

=cut

sub assign_fastq_files {
  my ( $self, $files ) = @_;
  return assign_files($files);
}

sub create_sorted_bam {
   my ( $self, $bam ) = @_;
   print "Creating sorted bam\n";
   my $cmd = $self->samtools . " sort ". $bam. " ";

   $bam =~ s/\.bam/\.sorted/;

   $cmd .= $bam;

   print $cmd,"\n";
   eval{
     `$cmd`;
   };

   if ($@) {
     throw( "Failed: $cmd" );      
   }


   my $out = $bam . "\.bam";
   print "Created $out\n";
   return $out;
}

sub collection_base_count{
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{collection_base_count} = $arg;
  }
  return $self->{collection_base_count};

}

sub skip_fragment {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{skip_fragment} = $arg;
  }
  return $self->{skip_fragment};
}

sub skip_mate_files {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{skip_mate_files} = $arg;
  }
  return $self->{skip_mate_files};
}
sub program_version {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{program_version} = $arg;
  }
  return $self->{program_version};
}


sub run_meta_info {

  my ( $self, $arg ) = @_;
  if ($arg) {

   if ( ! ($arg->isa("ReseqTrack::RunMetaInfo")) ){
     print Dumper ($arg), "\n";
     throw "Passed object other than RunMetaInfo\n"       	
   }
    $self->{run_meta_info } = $arg;
  }
  return $self->{run_meta_info };

}

sub percent_reads_used {
  my ( $self, $arg ) = @_;
  if ($arg) {
    $self->{percent_reads_used} = $arg;
  }
  return $self->{percent_reads_used}; 
}




1;

