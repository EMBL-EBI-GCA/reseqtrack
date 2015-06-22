#!/sw/arch/bin/perl -w
#RES cvs check

use warnings;
use strict;




use Getopt::Long;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::RunMetaInfoAdaptor;
use File::Basename;
use File::Path;
use Data::Dumper;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;


my $file;
my $file_list;
my $run        = 0;
my $verbose = 0;
my $mctr     = 0;
my $errors   = 0;
my $help    = 0;
my $expected_cols = 20;
my $count_cols = 1;
&GetOptions(
	    'dbhost=s'      => \$dbhost,
	    'dbname=s'      => \$dbname,
	    'dbuser=s'      => \$dbuser,
	    'dbpass=s'      => \$dbpass,
	    'dbport=s'      => \$dbport,
	    'file=s'            => \$file,
	    'file_list=s'      =>\$file_list,
	    'verbose'         => \$verbose,
	    'help'              => \$help,
	    'run'               => \$run,
	    'count_cols!'               => \$count_cols, 
	   );
if ($help) {
  useage();
}


my @bas_file;

#list of bas files
if ($file_list) {
  open (FH,'<',$file_list) || die "Failed to find $file_list";
  @bas_file = <FH>;
  close (FH);
}

if ($file) {
  push (@bas_file,$file);
}

die "No files specified" if (!@bas_file);



my $db = ReseqTrack::DBSQL::DBAdaptor->new(
					   -host   => $dbhost,
					   -user   => $dbuser,
					   -port   => $dbport,
					   -dbname => $dbname,
					   -pass   => $dbpass,
					  );


my $rmi_a = $db->get_RunMetaInfoAdaptor;




foreach my $inf (@bas_file) {
 
  chomp $inf;
  #  print $inf,"\n" ;
  $mctr   = 0;			# just a line counter
  $errors = 0;

  my $inf_name .= $inf;
  $inf_name = basename($inf_name);
  $inf_name =~ s/\.bam\.bas$//;
  # print "Col 1 name should be: $inf_name\n";



  open (FH,'<',$inf) || die "Failed to open: $file";

  while (<FH>) {

    next  if (/^bam/);

    my @data = split /\t/;

    my $cols_found =  scalar (@data);

     if ( ($cols_found != $expected_cols) && $count_cols){
       print "$mctr:Wrong number of columns: $cols_found should be 20\n"

    } 

    if ( $data[0] ne $inf_name) {
      print "$mctr:Name mismatch: $inf_name :  $data[0]\n";
    }



    my @sample_name = split /\./, $data[0]; #extract NAxxxxxx from col 1

    my  $meta_info = $rmi_a->fetch_by_run_id ( $data[6]);
    if ( ! -exists $meta_info->{run_id} ){
      print "$inf:Could not pull meta info for run_id = $data[6]\n";
      next;
    }

    $mctr++;

    if ($verbose) {
      print "$mctr\t";
      print "db [ " , $meta_info->{sample_name}, "\t",$meta_info->{run_id}, "\t",$meta_info->{library_name}." ]\t\t";
      print "bas file [ $sample_name[0]\t $data[6]\t   $data[5] ]\n";
	
    }
   
    #do sample names match?
    if ($meta_info->{sample_name} ne  $sample_name[0] ) {
      print "Error: ($mctr): sample  names bad ";
      print "for runid =    $meta_info->{run_id}::   bas file:  $sample_name[0]    run_meta_info: $meta_info->{sample_name}    $meta_info->{sample_name}\n";
      $errors++;
    }
       
    #do library names match ?
    if ($meta_info->{library_name} ne  $data[5] ) {
      print "Error: ($mctr): library names bad ";
      print "for runid =    $meta_info->{run_id}::   bas file:  $data[5]   run_meta_info: $meta_info->{library_name}    $meta_info->{sample_name}\n";
      $errors++;
    }
   

    foreach my $i (0 .. $#data) {
      if ( $data[$i] =~/unknown/) {
	print "Error 'unknown' entry in column $i: $data[$i]\n";
	$errors++;

      }
    }

    #    print Dumper ($meta_info);
  }
  close (FH);

  print "$inf Errors in this file = $errors\n" if ($errors);
  # print "$inf OK\n" if (! $errors);
     
}



sub useage{
  exec('perldoc', $0);
  exit(0);
}



=pod

=head1 NAME

ReseqTrack/scripts/qc/check_bas.pl

=head1 SYNOPSIS

This check to see if contents of bas file are consistent with run_meta_info
in g1k_archive_staging_track or other ReseqTrack db
Calls run_meta_info from database via "run_id" ( column 7 in bas file).
Checks sample names and library names are consistent with run_id specified in .bas file.

=head1 OPTIONS

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the
    standard port for mysql-g1kdcc.ebi.ac.uk
-file, input name of single bas file to check
-file_list, input list of bas files to check
-help, Binary flag to indicate if the help should be printed out
-verbose, watch line by line comparisons.


=head1 Examples

 perl $RESEQTRACK/scripts/check_bas.pl -file NA18525.454.MOSAIK.SRP000033.2009_11.bam.bas 
    -dbhost mysql-g1kdcc.ebi.ac.uk -dbport 4197 -dbuser g1krw -dbpass xxxxxxxxxx
    -dbname g1k_archive_staging_track

 perl $RESEQTRACK/scripts/check_bas.pl -file_list file.lst
    -dbhost mysql-g1kdcc.ebi.ac.uk -dbport 4197 -dbuser g1krw -dbpass xxxxxxxxx
    -dbname g1k_archive_staging_track

=cut




