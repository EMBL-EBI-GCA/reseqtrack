#!/sw/arch/bin/perl

use strict;
use warnings;

use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor; 
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::GeneralUtils;
use ReseqTrack::Tools::ArchiveUtils;
use ReseqTrack::Tools::AutoDiffTree;
use ReseqTrack::Tools::Loader::File;
use ReseqTrack::Tools::Loader::Archive;
use FileHandle;
use File::Path;
use File::Basename;
use Getopt::Long;

use Data::Dumper;

local $| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $host;


my $output_path;
my $verbose;
my $debug;
my $dir_to_tree = '/nfs/1000g-archive/vol1/ftp';
my $CHANGELOG = '/nfs/1000g-archive/vol1/ftp/CHANGELOG';
my $STAGING_DIR="/nfs/1000g-work/G1K/archive_staging/ftp";

my $check_old;
my $check_new;
my $old_tree_md5;
my $new_tree_md5;

my $old_tree_file;
my $new_tree_file;
my $log_dir;

my $test_mode;
my $skip_cleanup=0;
my $archive_sleep=240;
my $test=0;

&GetOptions( 
	    'dbhost=s'      => \$dbhost,
	    'dbname=s'      => \$dbname,
	    'dbuser=s'      => \$dbuser,
	    'dbpass=s'      => \$dbpass,
	    'dbport=s'      => \$dbport,
	    'host_name=s'        => \$host,
	    'dir_to_tree=s' => \$dir_to_tree,
	    'output_path=s' => \$output_path,	    
	    'verbose!'      => \$verbose, 
	    'check_old=s' =>\$check_old,
	    'check_new=s' =>\$check_new,
	    'staging_dir=s'   => \$STAGING_DIR,
	    'skip!'=>\$skip_cleanup,
	    'archive_sleep=s' => \$archive_sleep,
	    'debug!' => \$debug,
	    'test!' => \$test,
	    'log_dir=s' => \$log_dir,
	   );


throw("hostname not specified (-host). Required for Loader.pm") if (!$host);
mkpath($log_dir) unless(-d $log_dir);
throw("Can't run if ".$log_dir." log dir does not exist") unless(-d $log_dir);
my $logging_filepath = logging_filepath($log_dir);
my $log_fh = logging_fh($logging_filepath);

print $log_fh "Starting tree dump to:";
if ( !$output_path) {
  $new_tree_file = $STAGING_DIR . '/' . "current.tree";
  $new_tree_file =~ s/\/\//\//;
  print $log_fh $new_tree_file,"\n";	# if $verbose;
} else {
  $new_tree_file  = $output_path;
}





check_destinations($STAGING_DIR, $log_fh, $verbose);

my $dbA = ReseqTrack::DBSQL::DBAdaptor->new(
					    -host => $dbhost,
					    -user => $dbuser,
					    -port => $dbport,
					    -dbname => $dbname,
					    -pass => $dbpass,
					   );

throw "No DB connection established\n" if (! $dbA);


if ( ! ($check_old)  && ! ($check_new) ) {
  
  #fetch old tree object
  my $fa = $dbA->get_FileAdaptor;
  my $files = $fa->fetch_by_filename("current.tree");

  if (!$files || @$files == 0) {
    throw("Can't find the current.tree file in ".$dbname);
  } elsif (@$files >= 2) {
    throw("Seem to have ".@$files." versions of the current.tree file don't know ".
	  "what to do");
  }

  my $old_tree = $files->[0];

 

  throw("Old tree ".$old_tree->name." doesn't have a dbID in ".$dbname)
    unless($old_tree->dbID);

  $old_tree_md5 = $old_tree->md5;


  throw "Failed to get old_tree_md5\n" if (! $old_tree_md5 );


  
 
  dump_dirtree_summary($dir_to_tree,  $new_tree_file, undef, $fa);


  throw "Failed to get new_tree_md5\n" if (! -e  $new_tree_file );
  print $log_fh "Getting md5 of $new_tree_file\n";
  $new_tree_md5 = run_md5( $new_tree_file);

  
  $old_tree_file = $old_tree->name;

  print $log_fh "ftp current.tree md5     = ", $old_tree_md5,"\n";
  print $log_fh "Todays current.tree file = ", $new_tree_md5,"\n";

} else {
  # warning "Running in test mode, just comparing files\n";
  $test = 1;
  print $log_fh "check_new: $check_new\n";
  print $log_fh "check_old: $check_old\n";

  $new_tree_md5 = run_md5($check_new);
  $old_tree_md5 = run_md5($check_old);
  $new_tree_file      = $check_new;
  $old_tree_file      = $check_old;
  $test_mode = 1;
}





print $log_fh "new $new_tree_md5\n";
print $log_fh "old $old_tree_md5 \n";


if ($new_tree_md5 ne $old_tree_md5) {
  #create changelog files.
  my $DIFFTREE="ReseqTrack::Tools::AutoDiffTree";



  my $tree_diffs = $DIFFTREE->new(
				  -verbose       => $verbose,
				  -new_tree_file => $new_tree_file,
				  -old_tree_file => $old_tree_file,
				  -changelog     => $CHANGELOG,
				  -staging_dir   => $STAGING_DIR,
				 );

  $tree_diffs->create_log_files;
 


  my $files_to_process = $tree_diffs->files_to_archive_array();

 
  if (! $files_to_process ) {
    print $log_fh "No changelog files to process. current.tree should not have changed\n";
    unlink ( $new_tree_file);
    exit;
  }

  push (@$files_to_process , $new_tree_file );

  foreach my $file (@$files_to_process) {
    print $log_fh "process: $file\n";
    chmod (0755,$file);
  }  
 
  my $problem = 0;

  foreach my $file (@$files_to_process) {

    my $fa = $dbA->get_FileAdaptor;

    my $possible_existing = $fa->fetch_by_filename( basename($file) );

    if ( ($file =~ /CHANGELOG/) || ($file =~ /current\.tree/) ){

	 if (scalar (@$possible_existing) !=1){
	  print $log_fh "\nProblem:\n";
	  print $log_fh "$file should be in db.\n";
	  print $log_fh "Cannot find\n";
  	  $problem++;
 	 } 

       }
	else{
	 if (scalar (@$possible_existing) > 0){
  	  print $log_fh "Problem:\n";
	  print $log_fh "$file should not be in db.\n";
	  print $log_fh "It is dated changelog file and already exists\n";
          $problem++;
 	}
       } 

  }

  if ( $problem ){
    throw "Have got inconsistent file name problems output in ".$logging_filepath."\n";
  }


  if ( $test_mode) {
    print "In test mode comparing files. No database loading/archiving occurring\n";
    exit;

  }

 
  my $loader = ReseqTrack::Tools::Loader::File->new( 
						    -file      => $files_to_process,
						    -hostname  => $host,
						    -dbhost => $dbhost,
						    -dbname => $dbname,
						    -dbuser  => $dbuser,
						    -dbpass  => $dbpass,
						    -dbport  => $dbport,
						    -do_md5  => "1",
						    -update_existing=>"1",
						    -verbose=>$verbose,
						   );

  $loader->process_input();
  $loader->create_objects();
  $loader->sanity_check_objects();
  $loader->load_objects("1");



  my $archiver = ReseqTrack::Tools::Loader::Archive->new(
							 -db=>$dbA,
							 -file=> $files_to_process ,
							 -action => "archive",
							 -priority=>"90",
							 -max_number=>"1000",
							 -verbose=>$verbose,
							 -debug=>$debug, 
							);
  $archiver->process_input();
  $archiver->cleanup_archive_table($verbose);
  $archiver->sanity_check_objects();
  $archiver->archive_objects();


  my $max_tries = 10;
 
  if (!$skip_cleanup ) {
    print $log_fh "Starting final cleanup of archive table in $archive_sleep\n";
    sleep($archive_sleep);

    my $tries = 0;
    my $clean_archive_table = 1;

    while ($clean_archive_table) {
     
      my $obs_remaining = $archiver->cleanup_archive_table($verbose);
     
      if ($obs_remaining) {
	print $log_fh "Found $obs_remaining in archive table. Waiting $archive_sleep\n";
	sleep($archive_sleep);
	$tries++;
	
	if ($tries == $max_tries) {
	  print $log_fh "Tries $max_tries to clean up archive table. Giving up\n";
	  $clean_archive_table = 0;
	}
	
      } else {
	$clean_archive_table = 0;
      }
    }
  }
} else {
  print $log_fh "The tree files are exactly the same\n";
  unlink ($new_tree_file);
}




sub check_destinations {
  my ($staging_dir,$fh, $verbose) = @_;

  my $changelog_dir = $staging_dir . "/changelog_details";
  $changelog_dir =~ s/\/\//\//;

  print $fh $staging_dir,"\n"     if $verbose;
  print $fh $changelog_dir ,"\n"  if $verbose;

  throw "No staging_dir set" unless(-d $staging_dir);
  throw "$changelog_dir:\nAbove directory does not exist" unless(-d $changelog_dir);

  return;
}


sub logging_fh{
  my ($logging_filepath) = @_;
  my $open_string = ">".$logging_filepath;
  my $fh = FileHandle->new($open_string);
  unless($fh){
    throw("Failed to open ".$logging_filepath." $!");
  }
  return $fh;
}



sub logging_filepath{
  my ($log_dir) = @_;
  my $ident = int(rand(10000));
  my $date = current_date();
  my $tmp_name = "update_runmetainfo.".$date.".".$ident.".$$.log";
  return $log_dir."/".$tmp_name;
}

=pod

=head1 NAME

ReseqTrack/scripts/files/run_tree_for_ftp

=head1 SYNPOSIS

 This script should dump a 'current.tree' file into a staging directory
 then, compare it to the 'current.tree' file object in the database. If there are
 any differences the AutoDiffTree module should then compare the files
 and automatically construct changelog_details files and an amended CHANGELOG
 file in the staging directory. These files should then be automatically loaded
 into the specified database and archived.

=head1 OPTIONS

Database options

These set the parameters for the necessary database connection

 -dbhost, the name of the mysql-host
 -dbname, the name of the mysql database
 -dbuser, the name of the mysql user
 -dbpass, the database password if appropriate
 -dbport, the port the mysql instance is running on, this defaults to 4197 
          the standard port for mysql-g1kdcc.ebi.ac.uk

Standard options other than db paramters
 -host
 -dir_to_tree
 -staging_dir
 -output_path
 -skip
 -archive_sleep
 -verbose
 -debug

 -host            This is the name of the host which the filesystem 
                  is visible to so 1000genomes.ebi.ac.uk for ebi files.
                  Required by Loader module

 -dir_to_tree     directory to create current.tree file.
                  default: /nfs/1000g-work/G1K/archive_staging/ftp

 -staging_dir     directory location to create CHANGELOG and changelog_details
                  files. sub directory 'changelog_details' must exist.
                  default: /nfs/1000g-work/G1K/archive_staging/ftp

 -output_path    name for tree file to be generated.default 'current.tree'

 -skip           flag to skip archive table cleanup.

 -archive_sleep  wait period between each archive clean cycle.

 -test           create current.tree and changelogs but do not load/archive.
                 For testing purposes

 -verbose        usual what is going on output.
 -debug          more output


 Other test options
 For testing, you can just compare to tree files.
 -check_old     specifiy old tree file
 -check_new     specify  newer tree file

=head1 Examples

 $DB_OPTS= '-dbhost a_dbhost -dbport 4197 -dbuser a_user -dbpass ???? -dbname a_db -host a_host'
 
 Standard implementation
 perl  $Reseqtrack/scripts/file/run_tree_for_ftp.pl $DB_OPTS

 Test in tmp directory
 perl $Reseqtrack/scripts/file/run_tree_for_ftp.pl $DB_OPTS -staging_dir $PWD -host 1000genomes.ebi.ac.uk -test

 Compare tree files only
 perl $Reseqtrack/scripts/file/run_tree_for_ftp.pl $DB_OPTS  -check_old nov9_1212.tree -check_new $PWD/current.tree -staging_dir $PWD -host a_host 


=cut




