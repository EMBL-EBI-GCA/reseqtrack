#!/usr/local/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

use ReseqTrack::Tools::Exception qw(throw warning stack_trace_dump);
use ReseqTrack::Tools::FileUtils qw(create_objects_from_path_list create_history assign_type check_type);
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::IndexUtils qw(get_md5hash_from_sequence_index get_md5hash_from_alignment_index);
use ReseqTrack::DBSQL::DBAdaptor;

use ReseqTrack::File;
use ReseqTrack::Host;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $dir;
my @files;
my $list_file;
my $md5_file;
my $alignment_index_file;
my $sequence_index_file;
my $type;
my $remote;
my $host_name;
my $run;
my $verbose;
my $descend = 1;
my $die_for_problems = 1;
my $update_existing;
my $store_new;
my $use_dir;
my $check_md5;
my $md5_program = 'md5sum';
my $help;
my $assign_types = 1;
my $check_types = 1;
&GetOptions(
  'dbhost=s'       => \$dbhost,
  'dbname=s'       => \$dbname,
  'dbuser=s'       => \$dbuser,
  'dbpass=s'       => \$dbpass,
  'dbport=s'       => \$dbport,
  'dir=s' => \$dir,
  'use_dir!' => \$use_dir,
  'file=s@' => \@files,
  'list_file=s' => \$list_file,
  'md5_file=s' => \$md5_file,
  'alignment_index_file=s' => \$alignment_index_file,
  'sequence_index_file=s' => \$sequence_index_file,
  'type|file_type=s' => \$type,
  'host=s' => \$host_name,
  'run!' => \$run,
  'verbose!' => \$verbose,
  'descend!' => \$descend,
  'die_for_problems!' => \$die_for_problems,
  'update_existing!' => \$update_existing,
  'store_new!' => \$store_new,
  'check_md5!' => \$check_md5,
  'md5_program=s' => \$md5_program,
  'help!' => \$help,
  'assign_types!' => \$assign_types,
  'check_types!' => \$check_types,
    );

if($help){
  useage();
}

if(!$type && !$assign_types){
  throw("Must give load_files.pl a file type with -type");
}

if(!$host_name){
  throw("Must give load_files a host name with -host");
}

if(!$dbhost || !$dbname || !$dbuser){
  throw("Must provide database connection details with -dbhost -dbuser -dbpass ".
        "-dbport -dbname");
}

if($assign_types && !$type){
  $type = 'MUST_FIX';
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );

my $host = get_host_object($host_name, $db);

my @paths;

my $md5_hash;
if($dir){
  $dir =~ s/\/\//\//;
  $dir =~ s/\/$//;
}
if(@files){
  #file paths specified on commandline
  foreach my $file(@files){
    my $full_path = $file;
    if($use_dir && $dir){
      $full_path = $dir."/".$file;
      $full_path =~ s/\/\//\//g;
    }
    push(@paths, $full_path);
  }
}elsif($list_file){
  my $lines = get_lines_from_file($list_file);
  foreach my $line(@$lines){
    my $full_path = $line;
    if($use_dir && $dir){
      $full_path = $dir."/".$line;
      $full_path =~ s/\/\//\//g;
    }
    next if(-d $full_path);
    push(@paths, $full_path);
  }
}elsif($md5_file){
  $md5_hash = get_md5hash($md5_file);
}elsif($sequence_index_file){
  $md5_hash = get_md5hash_from_sequence_index($sequence_index_file);
}elsif($alignment_index_file){
  $md5_hash = get_md5hash_from_alignment_index($alignment_index_file);
}elsif($dir){
#get files from dir
  my ($paths, $hash) = list_files_in_dir($dir, 1);
  unless($descend){
    $paths = $hash->{$dir};
  }
  warning("No paths were found in ".$dir) if(!$paths || @$paths == 0);
  push(@paths, @$paths) if($paths && @$paths >= 1);
}else{
  throw("Must give the script options to describe what files to load ".
        "-file, -list_file, -md5_file, -index_file or -dir"); 
}

if($md5_hash && keys(%$md5_hash) >= 1){
  foreach my $key(keys(%$md5_hash)){
    my $md5 = $md5_hash->{$key};
    my $full_path = $key;
    if($use_dir && $dir){
      $full_path = $dir."/".$key;
      $full_path =~ s/\/\//\//g;
    }
    $md5_hash->{$full_path} = $md5;
    push(@paths, $full_path);
  }
}
print "Have ".@paths." file paths to input\n" if($verbose);
#create file objects
my $files = create_objects_from_path_list(\@paths, $type, $host);
$files = assign_type($files) if($assign_types);

if($check_types){
  my @wrong;
  foreach my $file(@$files){
    push(@wrong, $file) unless(check_type($file));
  }
  print STDERR "There are ".@wrong." files with the wrong type\n" if($verbose);
  foreach my $file(@wrong){
    print $file->name." ".$file->type." is wrong\n";
  }
  throw("There are problems with the file types") if(@wrong >= 1);
}
print "Have created ".@$files." files\n" if($verbose);

#sanity checks not sure if there should be any but check for existance if
#remote isn't specified
my @problems;
foreach my $file(@$files){
  my $string = $file->full_path." doesn't exist" unless($remote || 
                                                        -e $file->full_path);
  if(!$remote && -d $file->full_path){
    $string = $file->full_path." is a directory ";
  }
  push(@problems, $string) if($string);
}

if(@problems){
  foreach my $problem(@problems){
    print STDERR $problem."\n";
  }
  if($die_for_problems){
    throw(@problems." problems identified with input set dying");
  }
}
if($check_md5){
  warning("You are going to run md5s for ".@$files." files this may take a while")
      if(@$files >= 50);
}
my $fa = $db->get_FileAdaptor;

my @storage_problems;
 FILE:foreach my $file(@$files){
   if($check_md5){
     my $md5 = run_md5($file->full_path, $md5_program) if($run);
     $file->md5($md5);
   }
   if($md5_hash && keys(%$md5_hash)){
     my $md5 = $md5_hash->{$file->full_path};
     if($file->md5){
       unless($file->md5 eq $md5){
         print $file->full_path." has a different md5 to the once specified in ".
             "input Skipping the file\n";
         next FILE;
       }
     }
     $file->md5($md5);
   }
   eval{
     if($update_existing){
       my $existing = $fa->fetch_by_name($file->name);
       if($existing){
         $file->dbID($existing->dbID);
         my $history = create_history($file, $existing);
         $file->history($history) if($history);
         unless($history){
           next FILE;
         }
       }else{
         my $possible_existing = $fa->fetch_by_filename($file->filename);
         if($possible_existing){
           if(@$possible_existing == 1){
             my $existing = $possible_existing->[0];
             $file->dbID($existing->dbID);
             my $history = create_history($file, $existing);
             next FILE unless($history);
             $file->history($history) if($history);
           }elsif(@$possible_existing >= 2){
             my $for_update;
             foreach my $existing(@$possible_existing){
               if($existing->type eq $file->type){
                 if($for_update){
                   warning("Can't update ".$file->filename." there are multiple files ".
                           "which share its name and type");
                 }
                 $for_update = $existing;
               }
             }
             $file->dbID($for_update->dbID) if($for_update);
             my $history = create_history($file, $for_update) if($for_update);
             next FILE unless($history);
             $file->history($history) if($history);
             unless($for_update){
               print STDERR "There are ".@$possible_existing." possible existing files\n";
               foreach my $file(@$possible_existing){
                 print STDERR $file->dbID." ".$file->name."\n";
               }
               throw("Have multiple files linked to ".$file->filename." not sure how ".
                     "to update the file");
             }
           }
         }
       }
     }
     $fa->store($file, $update_existing, $store_new) if($run);
   };
   if($@){
     throw("Problem storing ".$file." ".$file->full_path." $@") if($die_for_problems);
     push(@storage_problems, $file->full_path." ".$@);
   }
}
foreach my $problem(@storage_problems){
  print STDERR $problem."\n";
}

sub useage{
  exec('perldoc', $0);
  exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/files/load_files.pl

=head1 SYNOPSIS

This script will take a list of file paths from a variety of options and load
the file table of a ReseqTrack database with the details

=head1 OPTIONS

Database options

This is the database the objects will be written to

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the standard
port for mysql-g1kdcc.ebi.ac.uk

File list options

-file, the path to a single file, this can appear on the commandline multiple times.

-list_file, a file containing a list of files paths, one to each line,

-md5_file, a file of output from md5sum in the format md5checkum filepath.

-alignment_index_file, a file in the format of the ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README.alignment.index, both the bam and the bai files will be stored

-sequence_index_file, a file in the format ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README.sequence.index.

-dir, If no other options are specified all the files are read from this directory
and the directory tree below it and stores. If you only want to store files in
that directory use -nodescend option.

-use_dir, If this is specified along with -dir string when any of the other file lists
are passed in the -dir string is associated with each path string given.

Other options

-type, This should be a string which will be associated with all files for example FASTQ if you are loading fastq files. There are no restrictions on what this is other 
than it should be shorter than 50characters, convention normally has each type in 
upper case and it is good if it is in someway informative about the files loaded. If you specify the -assign_types option you don't have to use this

-host, This is the name of the host which the filesystem is visible to so 1000genomes.ebi.ac.uk for ebi files.

-run, This is a binary option which needs no additional argument when specified the 
adaptor store method is run. By default this is off.

-verbose, This switchs on some print statements about the process. By default this is
off

-descend, If a directory is passed in this means the directory tree is descended and
all files below that point are stored. By default this option in on but you can use
the -nodescend option to prevent it happening

-die_for_problems, This is a binary option which alters how the script behaves when
it seems problems in the given file list. If this is switched on it prints out
all problems then throws and exception if not it just prints out the problems. By default this is off

-update_existing, This means if the store method finds a file with the same name it
will use the update method rather than the store method

-store_new, This means if the store method finds a file with the same name it will
still store a brand new line unless the md5s are the same

-check_md5, This makes the script run md5 checksums on each file

-md5_program, This is a string pointing to the md5sum program

-assign_types, This uses the ReseqTrack::Tools::FileUtils assign_type method to 
specify file types on the basis of the file path, this is on by default but it can be switched off by using -noassign_types

-check_types, This checks the type of each file object against an acceptable list, this is on by default but it can be switched off with -nocheck_types

-help, This makes the script print out its options

=head1 Examples

To load a single file

perl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname lec_test_track -file_type MISC -noassign_type -host 1000genomes.ebi.ac.uk -file /path/to/file -run

To load from an index file

perl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname lec_test_track  -host 1000genomes.ebi.ac.uk -sequence_index_file /path/to/sequence.index -run -use_dir -dir /root/to/sequence/files

To load from a directory

perl -dbhost mysql-g1kdcc -dbuser g1krw -dbpass XXXX -dbport 4197 -dbname lec_test_track -host 1000genomes.ebi.ac.uk -dir /path/to/dir -run


=cut
