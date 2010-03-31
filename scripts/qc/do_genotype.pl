use warnings;

use strict;


use Getopt::Long;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Host;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::FileAdaptor;
use ReseqTrack::Tools::BamUtils;
use FastQ;
use File::Basename;
use File::Path;
use File::stat;
use File::Copy;
use File::Temp qw/ tempfile tempdir /;
use Data::Dumper;



$| = 1;

my $GRC37  = "/nfs/1000g-work/G1K/work/reference/BWA/GRC37/human_g1k_v37";
my $g1k_work        = '/nfs/1000g-work/G1K/work';
my $samtools        =  '/nfs/1000g-work/G1K/work/bin/samtools/samtools';
my $bwa_exe         = '/nfs/1000g-work/G1K/work/bin/bwa-0.5.7/bwa';
my $bwamaleref      =  '/nfs/1000g-work/G1K/work/reference/BWA/NCBI36/human_b36_male';
my $bwafemaleref    =  '/nfs/1000g-work/G1K/work/reference/BWA/NCBI36/human_b36_female';

my $ncbimaleref     = '/nfs/1000g-work/G1K/work/reference/BWA/NCBI36/human_b36_male.fa';
my $ncbifemaleref   = '/nfs/1000g-work/G1K/work/reference/BWA/NCBI36/human_b36_female.fa';

my $glf             = "/nfs/1000g-work/G1K/work/bin/glftools/glfv3/glf";
my $snp_bin         = "/nfs/1000g-work/G1K/work/bin/ref/snps/hapmap3.snps.bin";

my $gender_file     = "/nfs/1000g-work/G1K/work/reference/BWA/NCBI36/sample_genders.csv";
my $snp_list        = "/homes/smithre/bin/glf_snps.lst";

my $staging_root = "/nobackup/vol1/smithre/GT_CHECK/STAGING/"; #processing directory 
my $internal_dir;		#process bam , create .bas in here
my $host_name;
my $file_type;
my $file;                       # unpaired fastq file
my $file1;			# assign_files = mate1
my $file2;			# assign_files = mate2
my $frag;
my $run_id;
my $type;
my $dbhost = "mysql-g1kdcc.ebi.ac.uk";
my $dbname ;
my $dbuser = "g1krw";
my $dbport = 4197;
my $dbpass = "thousandgenomes";
my $verbose;
my $run;
my $remote;
my $claimed_sample;
my $sample_label;
my $DEBUG = 0;
my $gender_hash;
my $ncbiref;
my $bwaref;

my @in_files;
my $subseq_files;
my $sub_size = 5000000; 
my @process_files;
my $grc37;


&GetOptions( 
	    'staging_root=s'   => \$staging_root,
	    'internal_dir=s'   => \$internal_dir,
	    'file1=s'          => \$file1,
	    'file2=s'          => \$file2, 
	    'frag=s'           => \$frag,
	    'run_id=s'         => \$run_id,
	    'type=s'           => \$type,
	    'host_name=s'      => \$host_name,
	    'run'              => \$run,
	    'claimed_sample=s' => \$claimed_sample,
	    'verbose'          => \$verbose,
	    'dbhost=s'         => \$dbhost,
	    'dbname=s'         => \$dbname,
	    'dbuser=s'         => \$dbuser,
	    'dbpass=s'         => \$dbpass,
	    'dbport=s'         => \$dbport,
	    'sub_size=s'       => \$sub_size,
	    'grc37'            => \$grc37,
	   );


$type="FILTERED_FASTQ";
$run=1;
$verbose =1;



my $G1K     = '/homes/smithre/OneKGenomes/modules/:';
my $RESEQ   = '/homes/smithre/OneKGenomes/ReseqTrack/modules/:';
my $CPAN    = '/homes/smithre/bin/perl/modules/:';

my $newperl5 =  $G1K . $RESEQ . $CPAN;





set_environ_vars ();


needed_files_check  ($samtools, $ncbimaleref, $ncbifemaleref,
		     $glf, $snp_bin, $gender_file,  $bwa_exe );


if (! $run) {
  print "\n\n** NOTE: Not in run mode. No important processing set to occur**\n\n";
}


#throw   "File names (-file1 -file2) not set"                           if (! $file1 && !$file2);
if (!$run_id) {

#  throw   "File name: $file1 does not appear to exist"          if (!-e  $file1);
#  throw   "File name: $file2 does not appear to exist"          if ($file2 && !-e  $file2);
}

throw   "Staging area location (-staging_dir) not set"         if (! $staging_root);
#throw   "Claimed Sample type does not specified"            if (! $claimed_sample);
#warning "Sub directory of processing dir not set (-internal_dir)" if (! $internal_dir);

my $rmi_a;

my $meta_info;
if ($run_id) {


  throw "File type not set" if (!$type);


  my $db = ReseqTrack::DBSQL::DBAdaptor->new(
					     -host => $dbhost,
					     -user => $dbuser,
					     -port => $dbport,
					     -dbname => $dbname,
					     -pass => $dbpass,
					    );



  # query information from Collection table
  my $ca         = $db->get_CollectionAdaptor;
  my $collection = $ca->fetch_by_name_and_type($run_id, $type);
  #print Dumper ($collection);
  if ( !defined ($collection)  ) {
    throw "Could not pull collection info for run_id = $run_id" ;
  }
 


  $rmi_a      = $db->get_RunMetaInfoAdaptor;
  $meta_info  = $rmi_a->fetch_by_run_id ($run_id);

  if (!defined ( $meta_info) ) {
    throw "Could not pull meta info for run_id = $run_id" ;
  }
   
 #  print Dumper ($collection);
 #  print Dumper ($meta_info);
 #  exit;

  print "run_meta_info  SAMPLE =  $meta_info->{sample_name}\n" if $verbose;
  $claimed_sample =  $meta_info->{sample_name};

  print  ref ($collection) if $verbose;

  throw("No collection object is found; please check the run_id and type\n") if (!$collection);

  my $others = $collection->others;
  #return a reference of an array of FileAdaptor objects.
  # A FileAdaptor object contains all info in a row in the File table

  throw("No File object is found for run_id $run_id\n") if (!$others);

  if ($verbose) { 
    my @coll_files;
    foreach my $q (@$others) {
      print $q->{name},"\t";
      print $q->{type},"\n";

      throw "Dead file $q->{name}" if ( ! -e  $q->{name} );
      push (@coll_files,  $q->{name});
    }

    print "@coll_files", "\n";
  }
  my  ($pe1 , $pe2, $frag) = assign_file_objects ($others);
  if ($pe1->{name} && $pe2->{name}) {
    print "$pe1->{name}     + $pe2->{name}\n" if $verbose;
  }
  print "$frag->{name}   by lonesome\n"   if ( $verbose && $frag->{name});


  $file1 = $pe1->{name}  if  $pe1->{name};
  $file2 = $pe2->{name}  if  $pe2->{name};
  $file  = $frag->{name} if  $frag->{name};

  if ( ($file1 && !$file2) || ($file2 && !$file1) ) {  
    print "Missing Paired end file";
    print "file 1 $file1\n";
    print "file 2 $file2\n";

  }

  if ($file1 && $file2) {
    my @tmp;
    push (@tmp, ($file1,$file2));
    push (@in_files, \@tmp);
  }


  if ( $file) {
    print "file $file\n" if $verbose;
    push (@in_files, [$file] );
  }

}


if (! $run_id) {		# giving file off command line
  my $db = ReseqTrack::DBSQL::DBAdaptor->new(
					     -host => $dbhost,
					     -user => $dbuser,
					     -port => $dbport,
					     -dbname => $dbname,
					     -pass => $dbpass,
					    );

  

  push (@in_files, $frag) if $frag; # unpaired fastq file

  push (@in_files, $file1) if  $file1;
  push (@in_files, $file2) if  $file2;
  throw "No files specified" if (! scalar @in_files);

}


my $snps = load_list_of_sample_snps ($snp_list); 

#print ref ($snps),"\n";

my @check = grep (/$claimed_sample/ ,@$snps);

print  "Failed   $run_id      $claimed_sample \n" if ( ! @check);
#print  "OK       $run_id      $claimed_sample    SNPS TO CHECK AGAINST\n";

#print @check,"\n";



# reference list of sample M/F,assign ncbi ref files
create_gender_hash ($gender_file,$meta_info->{instrument_platform} ); 

my $process_dir = '';

foreach my $r (@in_files) {
 
  my $now = time;

  foreach my $t ( @$r) {  
    print "Processing $t\n";
  }

  #where to create to genotype 
   $process_dir = create_process_dir ($staging_root, $run);
  print "Process directory =  $process_dir\n";


  chdir($process_dir) if $run;
 

  $subseq_files  = create_subsequence_sample (\@$r,$process_dir, $sub_size, $run);




  my $outsam     = from_fastq_to_sam   ($subseq_files,
					$bwaref, $process_dir,$run,
					$meta_info->{instrument_platform},$verbose);

  my $glf_output = run_glf_sample_check ( $outsam, $ncbiref,$glf ,$run , $verbose );

  my $passed     = check_glf_output ($glf_output, $claimed_sample, $verbose);

  $now = time - $now;
  printf("\n\nTotal running time: %02d:%02d:%02d\n\n",
	 int($now / 3600), int(($now % 3600) / 60),  int($now % 60));

  chdir($staging_root);
#  if ($passed) {
#    print "Passed genotype\n check:Deleting $process_dir\n";
  `chmod -R 775 $process_dir\*`;
  clean_up($process_dir);
#  }

  print "========================= Done =========================\n";
}
##### ## ## ##                         END            ## ##  ##  ######
#######################################################################
sub clean_up {

  my $pro_dir = shift;
  chdir($staging_root);
  `chmod -R 777 $pro_dir\*`;
  print "Deleting $pro_dir\n";
  rmdir($pro_dir);

}

sub  color_space_conversion{

 my $file = shift;
 my $IN;

 if ( $file =~ /.gz$/){
   print " XXXXXXXXXxxx zcating filepipe\n";
   open( $IN, "zcat $file |" ) or die "Can not open gz $file: $!\n";
 }

 if ( $file =~ /.fastq$/){
   open(  $IN, ,'<',"$file |" ) or die "Can not open fastq $file: $!\n";
 }

 throw ("No file handle") if (! $IN);

 (my $fh, my  $filename) = tempfile("g1k_XXXXX");

 print "\nTEMP FILE: $filename\n";

 while (<$IN>){
  
   my $line_id   = $_;
   my $line_seq  = <$IN>;
   my $line_sep  = <$IN>;
   my $line_qual = <$IN>;
   
   $line_seq = substr($line_seq, 2);
   $line_seq =~ tr/0123./ACGTN/;

   $line_qual = substr($line_qual, 2);

   if (length ($line_seq) != length ($line_qual)){
     die "Lengths do not match";
   }

   my $seq_len = length ($line_seq) -1 ;

   # $line_id =~ s/\d+$/$seq_len/;


   print $fh   $line_id, $line_seq, $line_sep , $line_qual ;



 }
  close ($fh);

 `gzip  $filename`;
 
 $filename .= ".gz";
 print "mv $filename $file\n";
 `mv $filename $file`;



}
 
sub create_subsequence_sample {
  my $infiles         = shift;
  my $process_dir     = shift;
  my $size            =shift;
  my $run             = shift if (@_);
  my @new_subsample_files;
  my $seq_list;

  foreach my $file (@$infiles) {

    my $outfile      = $$ . basename($file);

    print "sampling: #$file#\n" if $verbose;

    $seq_list = FastQ::sample( $file, $outfile,  $size)  if $run;
    #    FastQ::sample( $copy_fastq, $outfile, $size,  $seq_list)  if $run;


     if ($meta_info->{instrument_platform} =~ /.SOLID/i){
       print "HAVE TO DO SNEAKY CONSVERSION\n";
       color_space_conversion ($outfile);

     }




    push (@new_subsample_files, $outfile);

    if ($run) {
      my $tmp_size = stat($new_subsample_files[-1] )->size ;
      $tmp_size /= (1024*1000);
      printf "Size of $new_subsample_files[-1]  ~ %6.3f Mb (bytes)\n", $tmp_size ;
      print  "Warning: small sub sample size $tmp_size\n" if ($tmp_size < 1.0);
    }
  }
 # print "\n";

  my $ctr = 0;
  foreach my $file (@new_subsample_files) {
    $ctr++;
    print "$ctr   $file \n" if $verbose;
  }
  print "\n" if $verbose;
  return \@new_subsample_files;
}




sub check_glf_output {
  my $file   =shift;
  my $sample = shift;
  my $verbose = shift if (@_);
  my $top_score;
  my $top_sample;
  my %glf_scores;
  my %glf_rank;
  my $sec_score ;
  my $sec_sample;
  my $ctr = 0;
  my $i    = 0;
  my @all_scores; 

  print "Processing: $file\n" if $verbose;

  open (IN, '<', $file) || die "No glf output file to read";

  while (<IN>) {
    chomp;

    next if (m/^entropy/);	#skip first line

    my @line = split (/ /,$_);

    $line[1] =~ s/.snp// ;      #extract sample number

    my $samp   = $line[1];
    my $score  = $line[3];


   

    $all_scores [$i] = $line[3];
    $i++;

    if (! exists  $glf_scores{$samp}) {
      $glf_scores{$samp} = $score;
    } else {
      die "duplicate sample number in file";
    }

    if (! exists  $glf_rank{$samp}) {
      $glf_rank{$samp} = $ctr +1;
    } else {
      die "duplicate sample number in file";
    }


    $top_score = $score        if ($ctr==0);
    $top_sample = $line[1]     if ($ctr==0);

    $sec_score  = $score    if ($ctr==1);
    $sec_sample = $line[1]  if ($ctr==1);

    $ctr++;
  }

  if ( $top_score <= 0.0001) {
    print "Top score has likelihood of $top_score. Test probably failed\n"; 
    clean_up ($process_dir);
    warning "Test failed";
    return;
  }

  print "Top ranked hit = $top_sample ";
  print "    score  = $top_score\n";

  print "2nd ranked hit = $sec_sample ";
  print "    score  = $sec_score\n";


  my $ratio21 = $all_scores[1]/$all_scores[0];

  printf "Ratio of ranked 2 to 1              = %6.4f\n\n", 
    $ratio21;

  if (exists  $glf_scores{$sample}) {   

    my $ratio = $glf_scores{$sample}/$top_score;

    print "Claimed Sample:: $sample    Rank:: $glf_rank{$sample}\n";
   
    if ($glf_rank{$sample} ne "1"){
      print "\nClaimed Sample:: ";
      print "$glf_rank{$sample} :: $sample  $glf_scores{$sample}  Ratio = ";
      printf "%6.4f\n", $ratio;
    }
  }
  else{

    clean_up ($process_dir);
    throw "SAMPLE SNPS NOT THERE: $claimed_sample\n";
  }
  

  if (  $sample eq $top_sample ) {
    print "CLAIMED SAMPLE IS HIGHEST HIT\n";
  }


  my $sample_match = 0;
  if ($top_sample ne $sample) {
    print "Top sample does not match claimed sample\n\n";
  } else {
    $sample_match = 1;
  }

  if ( ($ratio21 > 1.2)   && $sample_match ) {
    print "\nConvincing match of data to Sample Number\n";
    return 1;
  }

  if ( ($ratio21 < 1.2) && ($ratio21 > 1.0)  && $sample_match ) {
    print "\nUnconvincing match of data to Sample Number\n" ;
    return 1;
  }

  
  if ( $ratio21 < 1.0 || !$sample_match ) {
    throw "WARNING: Sequence data appears not match Sample Number\n";
  }

}


sub run_glf_sample_check {

  my $insam   = shift;
  my $ncbiref = shift;
  my $glf     = shift;
  my $run     = shift;
  my $verbose = shift;
  print "sam file= $insam\n"   if $verbose;
  print "ncbiref  = $ncbiref\n" if $verbose;


  print "Creating bam:\n";
  print "samtools import $ncbiref  $insam $$.bam\n" if $verbose;
  my $make_bam = "samtools import $ncbiref  $insam $$.bam";
  eval {
    `$make_bam` if $run;
  };
  throw "$make_bam failed:$@" if ($@);

  # `samtools sort file.bam file.sorted` 
  #   samtools index file.sorted.bam     
  #   samtools pileup -g -f refseqs.fa file.sorted.bam > file.glf


  #**********************************************************

  #              Do bam QC right here

  #***********************************************************



  my $sort_bam=  "samtools sort  $$.bam  $$.sorted";
  print  "Running: $sort_bam\n" if $verbose;

  eval {
    `$sort_bam` if $run;
  };
  throw "$sort_bam failed:$@" if ($@);



 # print "\nsorted = $$.sorted\n";

  my $sorted  = $$ . ".sorted.bam";

  my $index_bam = " samtools index  $sorted ";
  print "\nRunning: $index_bam" if $verbose;

  eval {
    `$index_bam` if $run;
  };
 
  throw "$index_bam failed:$@" if ($@);



  my $make_glf="samtools pileup -g -f $ncbiref $$.sorted.bam > $$.glf"; 
  print    "\n$make_glf\n" if $verbose;

  eval {
    `$make_glf` if $run;
  };
  throw "$make_glf failed:$@" if ($@);

  my $glf_out = $$ . ".glf.out";

  my $run_glf =  $glf  . " checkGenotype $snp_bin $$.glf  > $glf_out";
  print $run_glf,"\n" if $verbose;

  eval {
    print "Running glf checkGenotype:$run_glf \n";
    `$run_glf` if $run;
  };
  throw "$run_glf failed:$@" if ($@);

  return $glf_out;

}




##############


sub from_fastq_to_sam {
  my $fastq_files     = shift;
  my $refseq          = shift;
  my $process_dir    = shift;
  my $run            = shift;
  my $platform       =shift;
  my $verbose         = shift if (@_);

  my @sai_files;
  my @create_sai_cmds;
  my $fastqgz      = '.fastq.gz';
  my $sai          = '.sai';
  my  $num_fq      = scalar @$fastq_files;
  my  $make_sam;


  foreach my $i (@$fastq_files) {
    print "processing: $i\n" if $verbose;
  }
  print "\n" if $verbose;

  #set bwa cmd line
  if ( $num_fq ==1) {		#single end
    $make_sam = "$bwa_exe samse  $refseq "; 
  }

  if ( $num_fq ==2) {		# paired ends
    $make_sam = "$bwa_exe sampe -a 2000 $refseq ";  
  }



  throw "\$make_sam not set. Have $num_fq fastq.gz files" if (! $make_sam);

  #create command to make .sai files from fastq files
  foreach my $i (@$fastq_files) {

    my $tmp_base = basename($i);
    my  $j = $i;
    $j =~ s/$fastqgz/$sai/;

    my $OPTS = '' ;

    $OPTS = " -cn  0.01 "   if (    $platform =~/SOLID/);
    $OPTS = " -q 20 -l 32 " if ( ! ($platform =~/SOLID/));


    push (@sai_files, $j);
    push (@create_sai_cmds, "$bwa_exe aln $OPTS $refseq $i 2>> ./log  >  $j ");
  }

  #create .sai files
  print "Creating sai files\n";
  my $ok=0;
  foreach my $i (@create_sai_cmds) {
    $ok = 0;
    print "$i\n" if $verbose;
    eval {
      `$i` if $run;
      $ok =  $? >> 8;
    };
    throw "$i failed:$ok" if ($ok);
  }

  #create cmd to produce .sam file
  foreach my $i (@sai_files) {
    $make_sam .=  "$i ";
  }


  foreach my $i (@$fastq_files) {
    $make_sam .=  "$i ";
  }
  my $outsam = $$ .".sam";
  $make_sam = $make_sam .  " 2> ./log  >  $outsam";


  print "\nCreating sam file \n";
  print $make_sam, "\n" if $verbose;

  eval {
    `$make_sam` if $run;
  };
  throw "$make_sam failed:$@" if ($@);

  throw "No .sam file ($outsam) created" if (!-e $outsam && $run); 
  throw "Failed to create sam file" if (  (-s  $outsam) ==0);
 
 return $outsam;
}

##########################

sub load_list_of_sample_snps {

  my $in =shift;

  open (my $IN ,"<", $in ) || die "Could not open snp list:$in";

  my @snps = <$IN>;

  close ($IN);

 # print "Total sample snps " , scalar (@snps),"\n";

  return \@snps;

}







sub create_gender_hash {
  my $ref_file = shift;
  my $platform = shift if (@_);
  my %genders;

  throw "Could not read gender reference list.\n$ref_file" if ( !-e $ref_file);

  open (IN, '<' ,$ref_file) || die "Failed to open $ref_file";

  while (<IN>) {
    chomp;
    my  @data = split /,/;
    if ( !$data[0]  || !$data[1] ) {
      throw "Incomplete gender designation";
    }

    if (!exists $genders{$data[0]}) {
      $genders{$data[0]} = $data[1];
    } else {
      warning "Duplicate gender entry";
    }

  }
  close (IN);


  # human_b36_male.fa or human_b36_female.fa
  if ( ! exists $genders{ $claimed_sample}) {
    throw ("Unknown sex for $claimed_sample\n"); 
    

  }
  else { 
   
    if ($genders{ $claimed_sample } eq "Male") {
      $ncbiref    = $ncbimaleref;
      $bwaref     = $bwamaleref; 
      print "Sample Male\n";
    }

    if ($genders{ $claimed_sample } eq "Female") { 
      $ncbiref    = $ncbifemaleref;
      $bwaref     = $bwafemaleref;
      print "Sample Female\n";
    }
  }

  print "ncbi ref = $ncbiref\n" if $verbose;
  print "bwa ref  = $bwaref\n"  if $verbose;

  print "PLATFORM IS $platform\n";
  if ($platform =~ /.SOLID/) {
    print "PLATFORM IS SOLID\n";

    my $cs = 'NCBI36/COLOR_SPACE';
    $bwaref  =~ s/NCBI36/$cs/;
   # $ncbiref =~ s/NCBI36/$cs/;




  }
  print "$bwaref\n";
  print "$ncbiref\n";

  return \%genders;
}   



sub create_tmp_bam_name {
  my $dir_name    = shift;
  my $file        = shift;
  my $verbose     = shift if (@_); 

  my $tmp_bam_name = $dir_name . '/' .  $$ . ".bam";
  print "Temporary bam file name=\n $tmp_bam_name\n" if $verbose;

  $tmp_bam_name=~ s/\/\//\//g;

  return  $tmp_bam_name;
}


sub needed_files_check{
  foreach my $f (@_) {
    print "-$f-\n" if ($DEBUG);
    throw "File not found: $f " if (!-e $f);
  }
}



sub set_environ_vars {

  #Needed to run samtools code
  if (! ($ENV{'PATH'} =~ /$g1k_work/) ) {

    $ENV{'PATH'}= 
      "$ENV{'PATH'}:$g1k_work";

    if (! ($ENV{'PATH'} =~ /$g1k_work/) ) {
      throw "samtools location not in \$PATH. Must append";

    }
  }
  print $ENV{'PATH'},"\n" if $DEBUG >1;

  if (! ($ENV{'PATH'} =~ /bwa/) ) {
    $ENV{'PATH'}= 
      "$ENV{'PATH'}:/nfs/1000g-work/G1K/work/bin/bwa-0.5.6";
    if (! ($ENV{'PATH'} =~ /bwa/) ) {
      throw "samtools location not in \$PATH. Must append";

    }
  }
}





sub  create_process_dir {
  my $root_dir  = shift;
  my $temp_dir  = "";
  my $run       = 0; 
  $run          = shift if (@_);


  return ""  if (! $run);

  my $tmp;
  $tmp = $root_dir . '/' ;
  $tmp =~ s/\/\//\//g;

  if ($tmp) {
    throw "Processing directory path looks incorrect.Need full path for safety"
      if (!( $tmp =~ /^\//) );  
  }


  $temp_dir = "NOT_RUNNING\n";
  $temp_dir =  tempdir ( DIR =>$tmp )  if $run;
 
  if (($run) && (! -e $temp_dir) ) {
    throw "$tmp: tmp directory does not exist and could not create: $@" ;
  }

  #print "Created  $temp_dir \n" if $run;
  `chmod 775 $temp_dir` if $run;
  return $temp_dir;
 
}
##





sub useage{
  exec('perldoc', $0);
  exit(0);
}



=pod

=head1 NAME

ReseqTrack/scripts/file/do_genotype.pl

=head1 SYNOPSIS

This script take a gzipped fastq file and sub-samples the contents. The sub-sample
is run through bwa and glftools as a genotype check to see if contents of fastq file
match the sample_name provided ( eg NA06984 is source of conents of fastq file

=head1 OPTIONS

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the
    standard port for mysql-g1kdcc.ebi.ac.uk
-help, Binary flag to indicate if the help should be printed out

-staging_root      #root directory to use for processing
-internal_dir       # sub-direcotry for processing in required ( If processing multple file for same sample).
-file1                  # First fastq file     (Labelled _1 for paired end files)
-file2                 # Second fastq file (Labelled _2 for paired end files)
-frag                 # unpaired fastq file
-verbose            # See what is happening
-run                   # Without you are basically in test mode
-claimed_sample;  # Sample name that the sequence data SHOULD originate  from
-sub_size          # set sub-sample size (default = 5000000 );

=head1 Examples


 NOTES: use FULL PATHS to fastq files
      : You have to specify claimed sample ( -claimed_sample NA06985) if specifying files
        on the command line. If using run_id the info is retreived from the DB.


perl /homes/smithre/OneKGenomes/reseq-personal/smithre/scripts/do_genotype.pl -staging_root $PWD  -frag $PWD/ERR000824.filt.fastq.gz  $G1K_DB_OPTIONS -claimed_sample NA11881-verbose -run -internal_dir fred

perl $ReseqTrack/scripts/qc/do_genotype.pl -staging_root $PWD  -file1 $PWD/ERR000824_1.filt.fastq.gz -file2 $PWD/ERR000824_2.filt.fastq.gz  $G1K_DB_OPTIONS -claimed_sample NA11881  -run -internal_dir fred

perl $ReseqTrack/scripts/qc/do_genotype.pl -staging_root $PWD  -file1 $PWD/ERR000824_1.filt.fastq.gz -file2 $PWD/ERR000824_2.filt.fastq.gz  $G1K_DB_OPTS -claimed_sample NA11881  -run -internal_dir fred


where $G1K_DB_OPTS = 
-dbhost mysql-g1kdcc.ebi.ac.uk -dbport 4197 -dbuser g1krw -dbpass xxxxx  -dbname g1k_archive_staging_track


=cut



