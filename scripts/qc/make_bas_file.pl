#!/sw/arch/bin/perl -w
use warnings;

use strict;

use Getopt::Long;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Host;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::FileAdaptor;
use ReseqTrack::Tools::BamUtils;


use File::Basename;
use File::Path;
use File::Copy;
use File::Temp qw/ tempfile tempdir /;
use Data::Dumper;


$| = 1;



my $NCBI36 = "/nfs/1000g-work/G1K/work/reference/BWA/NCBI36";
my $GRC37  = "/nfs/1000g-work/G1K/work/reference/BWA/GRC37";
my $maleref           = "/nfs/1000g-work/G1K/work/reference/BWA/NCBI36/human_b36_male.fa";
my $femaleref        = "/nfs/1000g-work/G1K/work/reference/BWA/NCBI36/human_b36_female.fa";
my $ref;
my $g1k_bin          = "/nfs/1000g-work/G1K/work/bin";
my $samtools_dir          = "/nfs/1000g-work/G1K/work/bin/samtools";
my $samtools_exe          = "/nfs/1000g-work/G1K/work/bin/samtools/samtools";
my $seqindex          = '/nfs/1000g-archive/vol1/ftp/pilot_data/pilot_data.sequence.index';
my $make_bas_file_cmd ='perl -MVertRes::Utils::Sam -e "VertRes::Utils::Sam->new->bas(shift, shift, shift)" ';
my $gender_file     = "/nfs/1000g-work/G1K/work/reference/sample_genders.csv";

my $G1K     = '/homes/smithre/OneKGenomes/modules/:';
my $RESEQ   = '/homes/smithre/OneKGenomes/ReseqTrack/modules/:';
my $CPAN    = '/homes/smithre/bin/perl/modules/:';
my $VertRes = '/nfs/1000g-work/G1K/work/bin/SangerVertRes/modules';
my $newperl5 =  $G1K .$RESEQ . $CPAN  . $VertRes;


my $ncbiref = "";

my $staging_root;      #directory where to put .bas file from cmd line
my $internal_dir;      #process bam , create .bas in here
my $final_type;  
my $host_name;
my $file_type;
my $file_name;
my $dbhost;
my $dbname;
my $dbuser;
my $dbport;
my $dbpass;
my $verbose;
my $run;
my $remote;
my $md5_program = 'md5sum';
my $have_tags;

#my $calc_md5=1;
#my $no_md5;
my $sample_label;
my $grc37;
my $sequenceindex_file;
my $release_date;
my $gender_hash;

my $store   = 0;
my $cleanup = 1;

my $DEBUG = 0;
$verbose =1 if ($DEBUG);


&GetOptions( 
	    'staging_root=s' => \$staging_root,
	    'internal_dir=s' => \$internal_dir,
	    'final_type=s'   => \$final_type,
	    'file_name=s'    => \$file_name,
	    'host_name=s'    => \$host_name,
	    'run'            => \$run,
	    'verbose'        => \$verbose,
	    'dbhost=s'       => \$dbhost,
	    'dbname=s'       => \$dbname,
	    'dbuser=s'       => \$dbuser,
	    'dbpass=s'       => \$dbpass,
	    'dbport=s'       => \$dbport,
	    'refseq=s'       => \$ncbiref,
	    'store'          => \$store,
	    'have_tags'   => \$have_tags,
	    'release_date=s' => \$release_date,
	    'grc37'          => \$grc37,
	    'si_file=s'      => \$sequenceindex_file,
#	   'no_md5'          => \$no_md5,    
	);



throw   "File name (-file_name) not set"                          if (! $file_name);
throw   "File name: $file_name does not appear to exist" if (!-e  $file_name);
throw   "Staging area location (-staging_dir) not set"            if (! $staging_root);
throw   "Final file type does not specified"                      if (! $final_type);
warning "Sub directory of processing dir not set (-internal_dir)" if (! $internal_dir);
die   "No release date" if ( ! $release_date);
die   "No sequence index file specified\n" if ( !$sequenceindex_file) ;


#$calc_md5 = 0 if ($no_md5); # if creating 


if (! $run){
  warning "-run option off. Will not generate .bas file";
}



#Needed to run Sanger code
if (! ($ENV{'PERL5LIB'} =~ /SangerVertRes/) ) {
  $ENV{'PERL5LIB'} = $newperl5;
  
  if (! ($ENV{'PERL5LIB'} =~ /SangerVertRes/) ) {
    throw "\$PATH does not include location of Sanger code\n";
  }
}
print $ENV{'PERL5LIB'},"\n" if $DEBUG > 1;


#Needed to run samtools code
#/nfs/1000g-work/G1K/work/bin/samtools
if (! ($ENV{'PATH'} =~ /$samtools_dir/) ) {

  $ENV{'PATH'}= 
    "$ENV{'PATH'}:$samtools_dir";

  if ( ! ($ENV{'PATH'} =~  /$samtools_dir/ )  ){
    throw "samtools location not in \$PATH. Must append";

  }
}
#print 'X' x 56, "\n";
print $ENV{'PATH'},"\n";# if $DEBUG >1;


    $ENV{'SAMTOOLS'}= $samtools_dir; 
    print   $ENV{'SAMTOOLS'},"\n";

print "\nProcessing: $file_name\n\n";

######################
#need md5 of original bam

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
   -host   => $dbhost,
   -user   => $dbuser,
   -port   => $dbport,
   -dbname => $dbname,
   -pass   => $dbpass,
  );

my $fa = $db->get_FileAdaptor();


my $md5;
my $g1k_entry;

$g1k_entry = $fa->fetch_all;

foreach my $h ( @{$g1k_entry}){
  next if ( ! $h->{md5} ) ;
  if ($h->{name} eq $file_name){
    print  $h->{name} , "  ", $h->{md5},"\n";
    $md5 = $h->{md5};
  }
}




#Check md5 of bam matches bam.md5 info.
if ( !$md5 ){
  print "\nmd5 not in db. Calculating md5\n";
  $md5 = run_md5($file_name, $md5_program) if ($run) ; 
}

if (!$md5){
  throw "Could not calculate md5 from $file_name" if $run;
}

$db->dbc->disconnect_when_inactive(1);

print "\nbam md5: ", $md5,"\n"  if $verbose;
#############


$sample_label =  extract_sample_label ($file_name, $verbose); # get NAxxxxx


$gender_hash = create_gender_hash ($gender_file); # reference list of sample M/F

if ($grc37){
  $ref = $GRC37;
}

# need to know if male/female for samtools reference fasta
# skip if using GRC37 build
if (! $ncbiref && !$grc37){


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



#where to create new .bas file ( should be in archive_staging somewhere)
my $process_dir = create_process_dir ($staging_root, $sample_label, $internal_dir,$run);
print "Process directory =  $process_dir\n" ;


my $file_base    = basename($file_name); #Target BAM file
print "File base name = $file_base\n\n"      if $DEBUG;


my $tmp_bam_name = create_tmp_bam_name ($process_dir, $file_base,$verbose); 


avoid_overwrite ( $tmp_bam_name,$verbose); #Do not overwrite existing bam files.


#process tmp bam file to produce .bas file
chdir($process_dir);
 my $ok=0;
 `pwd`;

my $tmp_bam ="$tmp_bam_name  2> /dev/null";
print  "$samtools_exe fillmd -b  $file_name $ref  2> /dev/null > $tmp_bam_name\n\n" if ($DEBUG || $verbose);


`$samtools_exe fillmd -b  $file_name $ref  2> /dev/null > $tmp_bam_name ` ;
$ok =  $? >> 8;
print " \$ok = $ok\n";
unlink ($tmp_bam_name) if ($ok);

throw "samtools failed: $@" if ($ok);


my $new_bas_file   = $file_base  . '.bas';
my  $run_make_bas = $make_bas_file_cmd . "$tmp_bam_name " . $sequenceindex_file .  "$process_dir/$new_bas_file";

print "Temporary bas file = $process_dir/$new_bas_file\n";
print  $run_make_bas,"\n\n" if ($DEBUG && $run);

eval{
  system ($run_make_bas)     if $run;
};
throw "VertRes failed: $@" if $@;


#put in correct ( original) bam md5 
correct_bas_file_convention ($md5 , $file_base, "$process_dir\/$new_bas_file") if $run;


=pod
if ( $store) {

  my $ha = $db->get_HostAdaptor;
  print "Create host object:$dbhost\n" if $verbose;
  my $host = $ha->fetch_by_name($dbhost);
  if (!$host) {
    $host = ReseqTrack::Host->new
      (
       -name => $dbhost,
       -remote => $remote
      );
  }




  if (!$db){
    my $db = ReseqTrack::DBSQL::DBAdaptor->new(
       -host   => $dbhost,
       -user   => $dbuser,
       -port   => $dbport,
       -dbname => $dbname,
       -pass   => $dbpass,
	      );
  }

  my $fax = $db->get_FileAdaptor;


  my $j = "$process_dir\/$new_bas_file";
  $j =~ s/\/\//\//g;


  print "Should store $j\n" if ($verbose || $DEBUG);



  my @file_array;
  push( @file_array, $j);

  my $files = create_objects_from_path_list(\@file_array , $final_type, $host);
  print "Have created ".@$files." file objects\n" if $verbose;



  foreach my $file (@$files) {
    print "Storing ", $file->full_path, "\n";
    $md5 = run_md5($file->full_path, $md5_program) ;
    print $md5,"\n";
    $file->md5($md5);
    warning "Store file in DB turned off";
    #  $fax->store($file) if $run;
  }

}


unlink($tmp_bam_name) if (-e $tmp_bam_name) ; 

if ($cleanup){
my $sanger_inline = $process_dir . '/'. "\_Inline";
$sanger_inline =~ s/\/\//\//g;

`chmod -R 775 $sanger_inline` if (-e $sanger_inline);
`rm -r $sanger_inline ` ;	#create by Sanger code.


if (-e "$file_name\.bas"){
  print "$file_name\.bas:\n exists could cross-compare\n";
}
}
=cut

print "$ file_name:    We are GOOD\n";

##### ## ## ##      END  ## ##  ##  ######
##########################################################################
##########################################################################
sub create_gender_hash {
  my $ref_file = shift;
  my %genders;

  throw "Could not read gender reference list.\n$ref_file" if ( !-e $ref_file);

  open (IN, '<' ,$ref_file) || die "Failed to open $ref_file";

  while (<IN>) {
    chomp;
   my  @data = split /,/;
    if  ( !$data[0]  || !$data[1] ){
      throw "Incomplete gender designation";
    }

    if (!exists $genders{$data[0]}){
      $genders{$data[0]} = $data[1];
    }
    else{
      warning "Duplicate gender entry";
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
  throw "Cannot not extract sample number from\n $file_name";
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
    throw "Processing directory path looks incorrect.Need full path for safety"
      if (!( $tmp =~ /^\//) );  
  }

  if (!-e $tmp) {
    print "Creating tmp dir in $tmp\n";
  `mkdir -p -m775 $tmp`;
    eval{
       mkdir ($tmp,0775) if $run;
     };

    throw "$tmp: tmp directory does not exist and could not create: $@" 
      if (! -e $tmp );
	print "Have created $tmp\n"; 
 }

  my $temp_dir =  tempdir ( DIR =>$tmp );
  if (! -e $temp_dir ){
    throw "$tmp: tmp directory does not exist and could not create: $@" ;
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


##
