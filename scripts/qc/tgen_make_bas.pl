use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use File::Path;
use File::Copy;
use File::Temp qw/ tempfile tempdir /;
use Data::Dumper;


$| = 1;

my $NCBI36 = "/nfs/1000g-work/G1K/work/reference/BWA/NCBI36";
my $GRC37  = "/nfs/1000g-work/G1K/work/reference/BWA/GRC37";
my $maleref           = $NCBI36."/human_b36_male.fa";
my $femaleref         = $NCBI36."/human_b36_female.fa";

my $samtools_dir      = "/nfs/1000g-work/G1K/work/bin/samtools";
my $samtools_exe      = "/nfs/1000g-work/G1K/work/bin/samtools/samtools";
my $seqindex          = '/nfs/1000g-archive/vol1/ftp/sequence_indices/20091216.sequence.index';

my $gender_file       = "/nfs/1000g-work/G1K/work/reference/sample_genders.csv";

my $have_md_tags      = 0;
my $make_bas_file_cmd = 'perl -MVertRes::Utils::Sam -e "VertRes::Utils::Sam->new->bas(shift, shift, shift)" ';

my $md5_program = 'md5sum';
my $G1K     = '/homes/smithre/OneKGenomes/modules/:';
my $RESEQ   = '/homes/smithre/OneKGenomes/ReseqTrack/modules/:';
my $CPAN    = '/homes/smithre/bin/perl/modules/:';
my $VertRes = '/nfs/1000g-work/G1K/work/bin/SangerVertRes/modules';
my $newperl5 =  $G1K .$RESEQ . $CPAN  . $VertRes;

my $ref;


my $ncbiref = "";


##

my $staging_root = ".";      #directory where to put .bas file from cmd line
my $internal_dir;      #process bam , create .bas in here
my $file_name;
my $verbose;
my $run;
my $grc37;

my $release_date;
my $md5;
my $gender_hash;
my $sample_label;
my $store   = 0;
my $cleanup = 1;
my $sequenceindex_file;


my $DEBUG = 1;
$verbose =1 if ($DEBUG);


&GetOptions( 
	    'staging_root=s' => \$staging_root,
	    'internal_dir=s' => \$internal_dir,
	    'file_name=s'    => \$file_name,
	    'run'            => \$run,
	    'verbose'        => \$verbose,
	    'refseq=s'       => \$ncbiref,
	    'md5=s'          => \$md5,
	    'sample_label=s' => \$sample_label,
	    'have_md_tags'   => \$have_md_tags,
	    'release_date=s' => \$release_date,
	    'grc37'          => \$grc37,
	    'si_file=s'      => \$sequenceindex_file;
	);



die   "File name (-file_name) not set"                          if (! $file_name);
die   "File name: $file_name does not appear to exist"          if (!-e  $file_name);
die   "Staging area location (-staging_dir) not set"            if (! $staging_root);
#die   "Sub directory of processing dir not set (-internal_dir)" if (! $internal_dir);
die   "No release date" if ( ! $release_date);
die   "No sequence index file specified\n" if ( !$sequenceindex_file) ;

if (! $run){
  print "-run option off. Will not generate .bas file";
}



#make sure required stuff is in $PATH
#set_envs ();



print "\nProcessing: $file_name\n\n";

######################
#need md5 of original bam, not tmp bam: goes in bas file
if ( !$md5 ){
  print "\nmd5 not in db. Calculating md5\n";
  $md5 = run_md5($file_name, $md5_program) if ($run) ; 
}

if (!$md5){
  die "Could n $? >> 8ot calculate md5 from $file_name" if $run;
}
print "\nbam md5: ", $md5,"\n"  if $verbose;
#############


#use to create tmp directory for processing.
if ( ! $sample_label){
  $sample_label =  extract_sample_label ($file_name, $verbose); # get NAxxxxx
}


if ($grc){
  $ref = $GRC37;
}


# need to know if male/female for samtools reference fasta
# skip if using GRC37 build
if (! $ncbiref && !$grc){


  $gender_hash = create_gender_hash ($gender_file); # reference list of sample M/F


  # human_b36_male.fa or human_b36_female.fa
  if (!exists $$gender_hash{ $sample_label }){
    print ("Unknown sex for $sample_label\n"); 

    my (%counts,$gender) = male_female_check($file_name,$verbose);
    $ref    = $maleref    if  ($gender eq "Male");
    $ref    = $femaleref  if  ($gender eq "Female");
  }
  else{ 
    $ref    = $maleref    if  ($$gender_hash{ $sample_label } eq "Male");
    $ref    = $femaleref  if  ($$gender_hash{ $sample_label } eq "Female");
  }
  print "Gender assigned\n" if ($ncbiref);

}

if (!$ref){
  die "Failed to assign reference to use ";
}

###########



#where to create tmp bam and  new .bas file. 
my $process_dir = create_process_dir ($staging_root, $sample_label, $internal_dir,$run);
print "Process directory =  $process_dir\n" ;




my $file_base    = basename($file_name); #Target BAM file
print "File base name = $file_base\n\n"      if $DEBUG;


my $tmp_bam_name = create_tmp_bam_name ($process_dir, $file_base,$verbose); 
avoid_overwrite ( $tmp_bam_name,$verbose); #Do not overwrite existing bam files.


#process tmp bam file to produce .bas file
chdir($process_dir);


my $ok=0;
print  "$samtools_exe fillmd -b  $file_name ref   > $tmp_bam_name\n\n" if ($DEBUG || $verbose);

`$samtools_exe fillmd -b  $file_name $ref  2> /dev/null  > $tmp_bam_name ` ;
$ok =  $? >> 8;
print " \$ok = $ok\n";
unlink ($tmp_bam_name) if ($ok);


print "\n\n\n CHECK $ok\n";

die "samtools failed: $@" if ($ok != 0);


my $new_bas_file   = $file_base  . '.bas';
my  $run_make_bas = $make_bas_file_cmd . "$tmp_bam_name " . $sequenceindex_file .  "$process_dir/$new_bas_file";



print "Temporary bas file = $process_dir/$new_bas_file\n";
print  $run_make_bas,"\n\n" if ($DEBUG && $run);

eval{
  system ($run_make_bas)     if $run;
};
die "VertRes failed: $@" if $@;


#put in correct ( original) bam md5 
correct_bas_file_convention ($md5 , $file_base, "$process_dir\/$new_bas_file") if $run;


print "$ file_name:    We are GOOD\n";

##### ## ## ##      END  ## ##  ##  ######
##########################################################################
##########################################################################
sub create_gender_hash {
  my $ref_file = shift;
  my %genders;

  die "Could not read gender reference list.\n$ref_file" if ( !-e $ref_file);

  open (IN, '<' ,$ref_file) || die "Failed to open $ref_file";

  while (<IN>) {
    chomp;
   my  @data = split /,/;
    if  ( !$data[0]  || !$data[1] ){
      die "Incomplete gender designation";
    }

    if (!exists $genders{$data[0]}){
      $genders{$data[0]} = $data[1];
    }
    else{
      die "Duplicate gender entry";
    }

  }
close (IN);

  return \%genders;
}   

sub extract_sample_label {
  my $file_name = shift;
  my $verbose     = shift if (@_); 
  print "\nextract sample number from :$file_name\n\n" if $verbose;

 #if given full path
 
 
    my @stuff =  split /[\/ | [.]/ , $file_name;
    foreach my $k (@stuff) {
      print " $k\n" if $verbose;
      if ( ($k =~ /[NA[\d][\d][\d][\d][\d]/)  ) {
	print "Sample Number:$k\n\n" if $verbose;
	return $k;
      }
    }
  die "Cannot not extract sample number from\n $file_name";
}
##





sub avoid_overwrite {
  #Do not overwrite existing bas files
  my $file_name = shift;
  my $verbose   = shift if (@_);

  if (-e "$file_name.bas" ) {
    print "$file_name.bas bas file already exists. Stopping\n"; 
    exit(0);
  }
  ;
  print "$file_name:\n Does not exist. Good to continue\n\n" if $DEBUG;
}
##



sub  create_process_dir {
  my $root_dir  = shift;
  my $sample   = shift;
  my $internal  = shift;
  my $run       = shift if (@_);

  my $tmp;
  $tmp = $root_dir . '/' . $sample. '/';
  $tmp .= $internal if $internal;
  $tmp =~ s/\/\//\//g;

  if ($tmp) {
    print $tmp,"\n";
    die "Processing dir path looks wrong.Need full path for safety"    
      if (!( $tmp =~ /^\//) );  
  }

  if (!-e $tmp) {
    print "Creating tmp dir in $tmp\n";
  `mkdir -p -m775 $tmp`;
    eval{
       mkdir ($tmp,0775) if $run;
     };

    die "$tmp: tmp directory does not exist and could not create: $@" 
      if (! -e $tmp );
	print "Have created $tmp\n"; 
 }

  my $temp_dir =  tempdir ( DIR =>$tmp );
  if (! -e $temp_dir ){
    die "$tmp: tmp directory does not exist and could not create: $@" ;
  }


  return $temp_dir;
}
##


sub correct_bas_file_convention{
  my $orig_bam_md5 = shift;
  my $bam_file = shift;
  my $bas_file  = shift;
  

  my @data;
  my @hold;
  my $newline;

  $bam_file =~   s/\.bam//g;

 # print ("fixing md5 with $orig_bam_md5\n");

  open (IN, '<' ,$bas_file) || die "Failed to open $bas_file";

  #Correct ill-formed bam file names
  while (<IN>) {
   # chomp;
   
    @data = split /\t/; 
  
    #ignore headers
    if ($data[0] eq "bam_filename"){
      $newline = join ("\t",@data);
        push (@hold, $newline);
      next;
    }


   
    if ( $data[0] =~ /^NA[\d][\d][\d][\d][\d]/) {
      print "correcting bam file name to  $bam_file\n";
      $data[0] = $bam_file;
    }

    $data[1] = $orig_bam_md5 ;
    
  #  if ( $data[4] =~ /ILLUMINA/) {
  #    $data[4] =  "SLX";
  #  }

    my $newline = join ( "\t", @data);

    push (@hold, $newline);
    
  }
  close (IN);



  my $tmp_file = $bas_file . ".org";
  move ($bas_file, $tmp_file);

  #Re-write file
  open (OUT, '>' ,"$bas_file") || die "Failed to open $bas_file for rewrite";
  foreach my $i (@hold) {
    print OUT $i;
  }
  close (OUT);
}
##


sub create_tmp_bam_name {
  my $dir_name    = shift;
  my $file             = shift;
  my $verbose      = shift if (@_); 

  my $tmp_bam_name = $dir_name . '/' .  $file. '_'. $$;
  print "Temporary bam file name=\n $tmp_bam_name\n\n" if $verbose;

  $tmp_bam_name=~ s/\/\//\//g;

  return  $tmp_bam_name;
}
##

sub male_female_check {
my $file    = shift;
my $verbose = shift if (@_);
my $sex;
my $ctr = 0;

die "$file:File does not exist" if (!-e $file);

if (! ($ENV{'PATH'} =~ /samtools/i) ) {
  die "samtools location not in \$PATH. Must append";
}

my @chr = ( 1..20, 'X', 'Y', "MT") ;


my %chrom_reads;
foreach (@chr){
  $chrom_reads{$_} = 0;
}

print "Checking designation of sample\n"          if $verbose;

open(BAMCAT,"samtools pileup $file |");
  my $i=0;
   while (<BAMCAT>) {
      chomp;
       my @data      = split(' ', $_);
       $chrom_reads { $data[0]}++;
   }
   close(BAMCAT);

if ($verbose){
  foreach my $key ( keys %chrom_reads){
    print  $key,"\t" ,$chrom_reads{$key},"\n"  if( $chrom_reads {$key} > 0) ;
  }
}

$sex = "male"    if ($chrom_reads{Y} > 0);
$sex = "female"  if ($chrom_reads{Y} == 0);


#print " in sub sex = $sex  $chrom_reads{Y}\n";

return  \%chrom_reads, $sex  ;

}##

sub set_envs{



}


sub run_md5{
  my ($file, $program) = @_;
  $program = 'md5sum' if(!$program);
  die("Can't run md5sum on a non existent file $file") unless (-e $file);
  my $cmd = $program." ".$file;
  open(FH, $cmd." | ") or die ("Failed to open ".$cmd);
  my $md5;
  while(<FH>){
    chomp;
    my @values = split;
    $md5 = $values[0];
  }
  close(FH);
  return $md5;
}


