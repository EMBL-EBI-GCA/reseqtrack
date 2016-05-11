#!/usr/bin/env perl

use strict;
use warnings;

use ReseqTrack::DBSQL::DBAdaptor; 
use ReseqTrack::Tools::CurrentTreeMaker;
use ReseqTrack::Tools::FileSystemUtils qw( get_lines_from_file);
use File::Basename qw(fileparse);
use Getopt::Long;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;

my $output_path;
my $dir_to_tree;
my $old_tree_name;
my $trim_dir;
my $file_list;
my %options;

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
      'trim_dir=s' => \$trim_dir,
      'options=s' => \%options,    
     );

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
              -host => $dbhost,
              -user => $dbuser,
              -port => $dbport,
              -dbname => $dbname,
              -pass => $dbpass,
             );

my $fa = $db->get_FileAdaptor;

if ($output_path && ! $old_tree_name)  {
  $old_tree_name = fileparse($output_path);
}

my $tree_maker = ReseqTrack::Tools::CurrentTreeMaker->new(
  -skip_regexes => $old_tree_name,
  -dir_to_tree => $dir_to_tree,
  -file_adaptor => $fa,
  -trim_dir => $trim_dir,
  -options => \%options,
);

if ($file_list) {
  my $files = get_lines_from_file($file_list);
  $tree_maker->files_to_tree($files);
}


$tree_maker->run;
if ($output_path) {
  $tree_maker->print_to_file($output_path);
}
else {
  $tree_maker->print_to_fh(*STDOUT);
}


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

 -trim_dir       Optional. Trims the full file path to a relative path
                 Default behaviour is to use dir_to_tree

 -file_list      Provide list of file in directory to dump tree for.
                 Side steps successive "find $input_dir" searches on large directories

 -old_tree_name  Optional. This file is not included in the output.  Default
                 is the basename of the output_path

 -verbose        usual what is going on output.
 -debug          more output



=head1 Examples

 $DB_OPTS= '-dbhost a_dbhost -dbport 4197 -dbuser a_user -dbpass ???? -dbname a_db -host a_host'
 
 Standard implementation
 perl  $Reseqtrack/scripts/file/dump_tree.pl $DB_OPTS \
  -dir_to_tree /nfs/1000g-archive/vol1/ftp/
  -output_path current.tree
  -old_tree_name current.tree

=cut




