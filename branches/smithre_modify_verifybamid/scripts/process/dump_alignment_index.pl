#!/sw/arch/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::DBSQL::DBAdaptor;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $type = 'BAM';
my $ftp_root = '/nfs/1000g-archive/vol1/ftp/';
my $output_file;
&GetOptions(
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'type=s' 		   => \$type,
  'ftp_root=s' 	   => \$ftp_root,
  'output_file=s'  => \$output_file,
    );


my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $ca = $db->get_CollectionAdaptor;

my $collections = $ca->fetch_by_type($type);
my $fa = $db->get_FileAdaptor;

my @sorted = sort{$a->name cmp $b->name} @$collections;
my $fh;
if($output_file){
  open(FH, ">".$output_file) or throw("Failed to open ".$output_file." $!");
  $fh = \*FH;
}else{
  $fh = \*STDERR;
}
#if ($type ne "EXOME_BI_BAM" && $type ne "EXOME_BCM_BAM" && $type ne "NCBI_BAM") { ### FIXME, change this when exome bi has bas files
	print $fh join("\t", ('BAM FILE', 'BAM MD5', 'BAI FILE', 'BAI MD5', 'BAS FILE', 'BAS MD5'))."\n" ;
#}
#else {
#	print $fh join("\t", ('BAM FILE', 'BAM MD5', 'BAI FILE', 'BAI MD5'))."\n" ;
#}
	
 COLLECTION:foreach my $collection(@sorted){
   #print $collection->name."\n";
   my @values = split /\./, $collection->name;
   my $sample_name = shift @values;
   my $analysis_grp = pop @values;
   my $pattern = $sample_name.".%.";
   $pattern .= join(".", @values) . "%" . $analysis_grp;
   
   my $files = $fa->fetch_all_like_name($pattern);
   my %hash;
   foreach my $file(@$files){
     if ($file->name !~ /$ftp_root/){
       print "Can't process ".$collection->name."  the file is not on the ftp\n";
       print $file->name."\n";
       #next COLLECTION;
       next;
     }
     
     if ($file->type =~ /INTERNAL_BA/) {
         print "Don't process internal files " . $file->name . "\n";
         next;
     }  
     
     if ($file->type =~ /WITHDRAWN/i) {
         print "Don't process withdrawn file " . $file->name . "\n";
         next;
     }
         
     next if ($file->name =~ /technical\/working\//); # won't include 
           
     #print $file->filename."\n";
     my $name = $file->filename;
     $name =~ s/\.bam//;
     $name =~ s/\.bai//;
     $name =~ s/\.bas//;
     push(@{$hash{$name}}, $file);
   }
   my @sorted_key = sort {$a cmp $b} keys(%hash);
   foreach my $key(@sorted_key){
    #print $key." ".scalar(@{$hash{$key}})."\n";
    my @files = @{$hash{$key}};
    my ($bam, $bas, $bai);
    foreach my $file(@files){
      if($file->filename =~ /\.bai$/){
        $bai = $file;
      }elsif($file->filename =~ /\.bas$/){
        $bas = $file;
      }elsif($file->filename =~ /\.bam$/){
        $bam = $file;
      }else{
        throw("Not sure what to do with ".$file->filename);
      }
    }
    if(!$bas || !$bai || !$bam){
    #if(!$bas && !$bai && !$bam){
      #throw("Not sure what to do for ".$key." don't have bas, bai or bam ".
      #      "files for it") if ($bam->type ne "EXOME_BI_BAM" && $bam->type ne "EXOME_BCM_BAM" && $bam->type ne "NCBI_BAM");  ### FIXME, change this when exome_bi_bam have bas files
      #throw("Not sure what to do for ".$key." don't have bas, bai or bam files for it"); 
      warning("Not sure what to do for ".$key." don't have bas, bai or bam files for it"); #FIXME
      next;
    }  
    #print "Printing index line for ".$key."\n";
    print_index_line($bam, $bas, $bai, $ftp_root, $fh);
   }
}


sub print_index_line{
  my ($bam, $bas, $bai, $ftp_root, $fh) = @_;
  my $bam_path = $bam->name;
  my $bam_md5 = $bam->md5;
  my $bai_path = $bai->name;
  my $bai_md5 = $bai->md5;
  my $bas_path = $bas->name;
  #my $bas_path = $bas->name if ($bam->type ne "EXOME_BI_BAM" && $bam->type ne "EXOME_BCM_BAM" && $bam->type ne "NCBI_BAM"); ### FIXME
  my $bas_md5 = $bas->md5;
  #my $bas_md5 = $bas->md5 if ($bam->type ne "EXOME_BI_BAM" && $bam->type ne "EXOME_BCM_BAM" && $bam->type ne "NCBI_BAM"); ### FIXME
  if($ftp_root){
    $bam_path =~ s/$ftp_root//;
    $bai_path =~ s/$ftp_root//;
    $bas_path =~ s/$ftp_root//;
  #  $bas_path =~ s/$ftp_root// if ($bam->type ne "EXOME_BI_BAM" && $bam->type ne "EXOME_BCM_BAM" && $bam->type ne "NCBI_BAM"); ### FIXME
  }
  #if ($bam->type ne "EXOME_BI_BAM" && $bam->type ne "EXOME_BCM_BAM" && $bam->type ne "NCBI_BAM") {### FIXME,
	  print $fh $bam_path."\t".$bam_md5."\t".$bai_path."\t".$bai_md5."\t".$bas_path.
    	  "\t".$bas_md5."\n";
  		#print $bam_path."\t".$bam_md5."\t".$bai_path."\t".$bai_md5."\t".$bas_path.
  		#    "\t".$bas_md5."\n";
  #}
  #else {
  #    print $fh $bam_path."\t".$bam_md5."\t".$bai_path."\t".$bai_md5."\n";
  #}    		
}

