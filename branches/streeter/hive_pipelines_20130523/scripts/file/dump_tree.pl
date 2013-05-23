#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::DBSQL::DBAdaptor; 
use ReseqTrack::Tools::FileSystemUtils qw( dump_dirtree_summary );
use Getopt::Long;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;

my $output_path;
my $dir_to_tree;
my $old_tree_name = 'current.tree';
my $output_prefix;
my $file_list;

&GetOptions( 
      'dbhost=s'      => \$dbhost,
      'dbname=s'      => \$dbname,
      'dbuser=s'      => \$dbuser,
      'dbpass=s'      => \$dbpass,
      'dbport=s'      => \$dbport,
      'dir_to_tree=s' => \$dir_to_tree,
      'output_path=s' => \$output_path,      
      'old_tree_name=s' => \$old_tree_name,      
      'file_list=s' => \$file_list,      
      'output_prefix=s' => \$output_prefix,      
     );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
              -host => $dbhost,
              -user => $dbuser,
              -port => $dbport,
              -dbname => $dbname,
              -pass => $dbpass,
             );

my $fa = $db->get_FileAdaptor;

if (!defined $output_prefix) {
  $dir_to_tree =~ /([^\/]+)\/*$/;
  $output_prefix = $1;
}

dump_dirtree_summary($dir_to_tree,  $output_path, $old_tree_name, $fa, $file_list, $output_prefix);


=pod

=head1 NAME

ReseqTrack/scripts/files/dump_tree.pl

=head1 SYNPOSIS

 This script creates a 'current.tree' file.  Nothing else.

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

 -dir_to_tree     directory to create current.tree file.

 -output_path    name for tree file to be generated. If not given, goes to stdout.

 -output_prefix  All files in the output_file will start with this prefix.
                 Default is the name of dir_to_tree e.g. 'ftp' (not the full path)

 -file_list      Provide list of file in directory to dump tree for.
                 Side steps successive "find $input_dir" searches on large directories

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




