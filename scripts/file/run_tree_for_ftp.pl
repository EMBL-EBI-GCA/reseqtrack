#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::Tools::Exception qw();
use ReseqTrack::DBSQL::DBAdaptor; 
use ReseqTrack::Tools::FileSystemUtils qw(check_file_exists check_directory_exists);
use ReseqTrack::Tools::GeneralUtils qw(current_date current_time);
use ReseqTrack::Tools::Loader::File;
use ReseqTrack::Tools::Loader::Archive;
use ReseqTrack::Tools::CurrentTreeMaker;
use ReseqTrack::Tools::CurrentTreeDiffer;
use File::Basename qw(fileparse);
use Getopt::Long;

local $| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $host = '1000genomes.ebi.ac.uk';


my $changelog_name = 'CHANGELOG';
my $changelog_details_name = 'changelog_details';
my $file_tree_name = 'current.tree';
my $dir_to_tree = '/nfs/1000g-archive/vol1/ftp';
my $staging_dir="/nfs/1000g-work/G1K/archive_staging/ftp";
my %options;
my $old_tree_dir;
my $old_changelog_dir;
my $log_dir;
my $trim_dir;

my $skip_cleanup=0;
my $skip_archive=0;
my $skip_load=0;
my $archive_sleep=240;

my $command_line = join(' ', $0, @ARGV);
&GetOptions( 
        'dbhost=s'      => \$dbhost,
        'dbname=s'      => \$dbname,
        'dbuser=s'      => \$dbuser,
        'dbpass=s'      => \$dbpass,
        'dbport=s'      => \$dbport,
        'host_name=s'        => \$host,
        'dir_to_tree=s' => \$dir_to_tree,
        'staging_dir=s'   => \$staging_dir,
        'old_tree_dir=s' => \$old_tree_dir,
        'file_tree_name=s'   => \$file_tree_name,
        'old_changelog_dir=s' => \$old_changelog_dir,
        'changelog_name=s'   => \$changelog_name,
        'skip!'=>\$skip_cleanup,
        'skip_load!'=>\$skip_load,
        'skip_archive!'=>\$skip_archive,
        'archive_sleep=s' => \$archive_sleep,
        'log_dir=s' => \$log_dir,
        'options=s' => \%options,
       );
$old_tree_dir //= $dir_to_tree;
$old_changelog_dir //= $dir_to_tree;

my $old_tree_path = $old_tree_dir . '/' . $file_tree_name;
my $new_tree_path = $staging_dir . '/' . $file_tree_name;
$old_tree_path =~ s{//}{/}g;
$new_tree_path =~ s{//}{/}g;
check_file_exists($old_tree_path);


my $log_fh;
if ($log_dir) {
  check_directory_exists($log_dir);
  my $date = current_date;
  $log_fh = File::Temp->new( UNLINK => 0, DIR => $log_dir,
    TEMPLATE => "run_tree_for_ftp.$date.XXXX", SUFFIX => '.log');
  #print "Writing log messages to ", $log_fh->filename, "\n";
  print_log($log_fh, $command_line, "\n");
}

# create the tree diffs object early in the script to make sure date is correct
# because the cron job is run close to midnight.
my $tree_diffs = ReseqTrack::Tools::CurrentTreeDiffer->new(
                -old_tree => $old_tree_path,
                -old_changelog_dir => $old_changelog_dir,
                -changelog_name     => $changelog_name,
                -changelog_details_name => $changelog_details_name,
                -output_dir   => $staging_dir,
                -skip_regexes => [qr{$changelog_name$}, qr{/$changelog_details_name/}],
               );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
                        -host => $dbhost,
                        -user => $dbuser,
                        -port => $dbport,
                        -dbname => $dbname,
                        -pass => $dbpass,
                       );
throw("No DB connection established", $log_fh) if (! $db);

my $clean_archiver  = ReseqTrack::Tools::Loader::Archive->new(
                                  -db=>$db,
                                  -no_lock=>1,
                                  -action=>'archive',
                                 );

print_log($log_fh, "checking archive table is clean\n");
if ($clean_archiver->cleanup_archive_table()) {
  throw("Stopping tree process as there are still objects in the archive table and this may cause issues", $log_fh);
}

my $fa = $db->get_FileAdaptor;

my $tree_maker = ReseqTrack::Tools::CurrentTreeMaker->new(
  -skip_regexes => $file_tree_name,
  -dir_to_tree => $dir_to_tree,
  -file_adaptor => $fa,
  -trim_dir => $trim_dir,
  -options => \%options,
);

$db->dbc->disconnect_when_inactive(2);

print_log($log_fh, "starting to build file tree\n");
$tree_maker->run;
print_log($log_fh, "dumping file tree to $new_tree_path\n");
$tree_maker->print_to_file($new_tree_path);
$tree_maker->print_errors_to_fh();
$tree_maker->print_errors_to_fh($log_fh) if $log_fh;

$tree_diffs->new_tree($new_tree_path);
if (! $tree_diffs->quick_diff()) {
  print_log($log_fh, "Old and new tree files are identical.  Nothing more to do\n");
  unlink($new_tree_path);
  exit(0);
}

#create changelog files.

my $ftra = $db->get_FileTypeRuleAdaptor;
my $file_type_rules = $ftra->fetch_all_in_order;

$tree_diffs->file_type_rules($file_type_rules);

print_log($log_fh, "starting to build file changelogs\n");
$tree_diffs->run;
print_log($log_fh, "dumping changelogs\n");
$tree_diffs->print_changelogs;
$tree_diffs->print_errors_to_fh();
$tree_diffs->print_errors_to_fh($log_fh) if $log_fh;
my $new_changelog_files = $tree_diffs->output_files;
print_log($log_fh, "number of new changelog files: ", scalar @$new_changelog_files, "\n");

if ($skip_load) {
  print_log($log_fh, "Not loading files.  Nothing more to do\n");
  exit(0);
}
print_log($log_fh, join("\n", "Loading files:", $new_tree_path, @$new_changelog_files), "\n");

FILE:
foreach my $file (@$new_changelog_files) {
  my $basename = fileparse($file);
  next FILE if $basename eq $changelog_name;
  my $existing = $fa->fetch_by_filename($basename);
  if (@$existing) {
    throw("File with name $basename already exists in the database: " .  join(',', map {$_->name} @$existing), $log_fh);
  }
}



my $loader = ReseqTrack::Tools::Loader::File->new( 
                          -file      => [@$new_changelog_files, $new_tree_path],
                          -hostname  => $host,
                          -dbhost => $dbhost,
                          -dbname => $dbname,
                          -dbuser  => $dbuser,
                          -dbpass  => $dbpass,
                          -dbport  => $dbport,
                          -do_md5  => "1",
                          -update_existing=>"1",
                         );

$loader->process_input();
$loader->create_objects();
$loader->sanity_check_objects();
$loader->load_objects("1");

if ($skip_archive) {
  print_log($log_fh, "Not archiving files.  Nothing more to do\n");
  exit(0);
}
print_log($log_fh, join("\n", "Archiving files:", $new_tree_path, @$new_changelog_files), "\n");

my $archiver = ReseqTrack::Tools::Loader::Archive->new(
                           -db=>$db,
                           -file => [@$new_changelog_files, $new_tree_path],
                           -action => "archive",
                           -priority=>"90",
                           -max_number=>"1000",
                          );
$archiver->process_input();
$archiver->cleanup_archive_table();
$archiver->sanity_check_objects();
$archiver->archive_objects();


if ($skip_cleanup) {
  print_log($log_fh,"Not waiting to cleanup archive table.  Nothing more to do\n");
  exit(0);
}

print_log($log_fh, "Starting final cleanup of archive table in $archive_sleep\n");

my $max_tries = 10;
my $tries = 0;
CLEAN:
while (1) {
  $tries += 1;
  sleep($archive_sleep);
  my $obs_remaining = $archiver->cleanup_archive_table();
  last CLEAN if !$obs_remaining;
  next CLEAN if $tries < $max_tries;
  print_log($log_fh, "Tries $tries to clean up archive table. Giving up\n");
  last CLEAN;
}


sub throw {
  my ($msg, $log_fh) = @_;
  print_log($log_fh, $msg, "\n");
  close($log_fh);
  ReseqTrack::Tools::Exception::throw($msg);
}

sub print_log {
  my ($log_fh, @args) = @_;
  return if ! $log_fh;
  print $log_fh current_time(), "\t", @args;
  $log_fh->flush;
}


=pod

=head1 NAME

ReseqTrack/scripts/files/run_tree_for_ftp

=head1 SYNPOSIS

 This script should dump a 'current.tree' file into a staging directory
 then, compare it to the 'current.tree' file object in the database. If there are
 any differences the CurrentTreeDiffer module should then compare the files
 and automatically construct changelog_details files and an amended CHANGELOG
 file in the staging directory. These files should then be automatically loaded
 into the specified database and archived.

 Overriding default behaviour:
 The modules CurrentTreeMaker and CurrentTreeDiffer have been written so that
 other projects can override the default behaviour.

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

 -host            Default is 1000genomes.ebi.ac.uk
                  Required by Loader module

 -options         for constructing a hash of options for the CurrentTreeMaker module
                  e.g. -options skip_base_directories=1

 -dir_to_tree     directory to create current.tree file.
                  default: /nfs/1000g-work/G1K/archive_staging/ftp

 -staging_dir     Files are written out to this directory before archiving
                  default: /nfs/1000g-work/G1K/archive_staging/ftp

 -old_tree_dir    Default is to be equal to dir_to_tree
                  The old tree file must be present in this directory

 -file_tree_name  Default is current.tree
                  The output file has this name
                  This file must also be present in the old_tree_dir

 -old_changelog_dir Default is to be equal to dir_to_tree
                  The old changelog file must be present in this directory

 -changelog_name  Default is CHANGELOG
                  The output file has this name
                  This file must also be present in the old_tree_dir

 -skip           Default 0. Flag to skip archive table cleanup.

 -skip_load      Default 0. Flag to skip loading of output files into the database

 -skip_archive   Default 0. Flag to archiving of files

 -archive_sleep  wait period between each archive clean cycle.

 -log_dir        Optional. Log messages are written to this directory


=head1 Examples

 $DB_OPTS= '-dbhost a_dbhost -dbport 4197 -dbuser a_user -dbpass ???? -dbname a_db -host a_host'
 
 Run it like this for the thousandgenomes project:
 perl  $Reseqtrack/scripts/file/run_tree_for_ftp.pl $DB_OPTS -options skip_base_directories=1

 For any other projecct:
 perl  $Reseqtrack/scripts/file/run_tree_for_ftp.pl $DB_OPTS \
        -dir_to_tree /path/to/archive -staging_dir /path/to/staging


=cut




