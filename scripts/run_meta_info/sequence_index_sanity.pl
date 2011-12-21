#!/sw/arch/bin/perl -w

use strict;

use Getopt::Long;

use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Tools::SequenceIndexSanity;
use ReseqTrack::Tools::RunMetaInfoUtils;

$| = 1;

my $index_file;
my $ftp_root;
my $syntax;
my $filesystem;
my $verbose;
my $file_exists;
my $ftp_path_must_contain = 'sequence_read';
my $dir_must_contain = 'ftp/data';
my $help;
my $log;
my @command_args = @ARGV;


#some default time for logging purpose:

    my($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = gmtime();
    my $year = 1900 + $yearOffset;
    my $theGMTime = "$year-$dayOfMonth-$month $hour:$minute:$second";
    #example: 2009-3-20, 16:12:43

    

&GetOptions( 
  'index_file=s' => \$index_file,
  'ftp_root=s' => \$ftp_root,
  'verbose!' => \$verbose,
  'syntax!' => \$syntax,
  'filesystem!' => \$filesystem,
  'log' => \$log,
  'ftp_path_must_contain=s' => \$ftp_path_must_contain,
   'dir_must_contain=s' => \$dir_must_contain,
  'help!' => \$help,
    )or useage(@command_args);


if(!$syntax && !$filesystem){
  warning("Can't run any sanity checks without specifying at least one of ".
        "-syntax or -filesystem");
  $help = 1;
}

if(!$index_file || !(-e $index_file)){
  warning("Can't sanity check an index file ".$index_file." which doesn't exist");
  $help = 1;
}


if($help){
  useage(@command_args);
}

my $index_hash = get_index_hash($index_file);


if($filesystem){
	print "Filesystem check:\n";
  	print "Comparing ".$index_file." to ".$ftp_root."\n";
  	my ($not_index, $not_dir) = compare_index_to_ftp($index_file, $ftp_root, 
                                                   $ftp_path_must_contain, $dir_must_contain);
  	print "There are ".@{$not_index->list}." files not in the index.\n";
  	print "There are ".@{$not_dir->list}." files not on the ftp.\n";

  	my ($FTP_total_count,$FTP_faulty_count,%FTP_log) = check_index_file_existance_to_ftp($index_file, $ftp_root);
 	print "FTP Check:\t".$FTP_faulty_count." are faulty of a total of ".$FTP_total_count." rows.\n";
	
	my ($PAIRED_total_count,$PAIRED_faulty_count,%PAIRED_log) = check_paired_files($index_file, $ftp_root);
  	print "PAIRED Check:\t".$PAIRED_faulty_count." files (up to 2 per row) are faulty of a total of ".$PAIRED_total_count." rows.\n";
  
  	
  	if($verbose){
  		print "Files not on index\n";
    	foreach my $file(@{$not_index->list}){
    		print $file."\n";
    	}
    	print "\nFiles not on filesystem ".$ftp_root." according to intersection ".
            "check \n";
    	foreach my $file(@{$not_dir->list}){
    		print $file."\n";
    	}
    	print "\nFiles not on FTP\n";
   		foreach my $file(keys(%FTP_log))
  		{
  			print $file."\n";
  		}
  	
  		print "\nFiles not correctly paired:\n";
   		foreach my $file(keys(%PAIRED_log))
  		{
  			print $file."\n";
  		}


  	}
  	
  	#START LOG FILE WRITING
	open(DAT,">>log_filesystem.csv") || die("Cannot Open File");
	print DAT "$theGMTime, $FTP_total_count, $FTP_faulty_count, $PAIRED_total_count, $PAIRED_faulty_count,".@{$not_index->list}.",".@{$not_dir->list}."\n";
	close(DAT);
	if($log)
	{
		
		#DETAILED LOG FILES
		#Contain Error \tab\ filename

		#FTP error log:
		if($FTP_faulty_count > 0)
		{
			#error log
			open(DAT, ">log_filesystem_FTP_$theGMTime.txt");
			print DAT "ERROR\t File\n";
			foreach my $file(keys(%FTP_log))
	  		{
	  			print DAT $file."\n";
	  		}
			close(DAT);
		}
		
		#PAIRED error log:
		if($PAIRED_faulty_count > 0)
		{
			#error log
			open(DAT, ">log_filesystem_PAIRED_$theGMTime.txt");
			print DAT "ERROR\t File\n";
			foreach my $file(keys(%PAIRED_log))
	  		{
	  			print DAT $file."\n";
	  		}
			close(DAT);
		}
		
		#Missing on INDEX error log:
		if(@{$not_index->list} > 0)
		{
			#error log
			open(DAT, ">log_filesystem_INDEX_$theGMTime.txt");
			print DAT "ERROR\t File\n";
    		foreach my $file(@{$not_index->list}){
    			print DAT "Missing file:\t".$file."\n";
    		}
			close(DAT);
		}
		
		#Missing on FTP error log:
		if(@{$not_dir->list} > 0)
		{
			#error log
			open(DAT, ">log_filesystem_FTP2_$theGMTime.txt");
			print DAT "ERROR\t File\n";
    		foreach my $file(@{$not_dir->list}){
    			print DAT "Missing file:\t".$file."\n";
    		}
			close(DAT);
		}
	}
}


if($syntax){

  print "Syntax check:\n";
  check_header($index_file);
  my $column_count_hash = check_column_count_sanity($index_file);
  print "There are ".keys(%$column_count_hash)." lines in ".$index_file.
      " with the wrong number of columns\n" if($column_count_hash && 
                                               keys(%$column_count_hash) >= 1);
  if($verbose){
    my $count;
    foreach my $key(keys(%$column_count_hash)){ 
      print $key." ".$column_count_hash->{$key}."\n";
      $count++;
      last if($count >= 10);
    }
  }
  my $column_value_hash = has_empty_columns($index_file);
  print "There are ".keys(%$column_value_hash)." lines in ".$index_file.
      " with undefined columns\n" 
      if($column_value_hash && keys(%$column_value_hash) >= 1);
  if($verbose){
    if($column_value_hash && keys(%$column_value_hash)){
      my $count;
    KEY:foreach my $key(keys(%$column_value_hash)){
      my $values = $column_value_hash->{$key};
      foreach my $index(@$values){
        print "We have a problem for ".$key." as ".$index." is undefined\n"
            if($index == 0 || $index == 1 || $index == 2 || $index == 3 ||
               $index == 5 || $index == 6 || $index == 10 || $index == 14 ||
               $index == 18);
        $count++;
        last KEY if($count >= 10);
      }
    }
    }
  }
  my $has_nulls_hash = has_null_strings($index_file);
  print "There are ".keys(%$has_nulls_hash)." lines in ".$index_file." with ".
      "NULL strings in columns\n";
  if($verbose){
    my $count;
    foreach my $key(keys(%$has_nulls_hash)){
      my $entries = $has_nulls_hash->{$key};
      print $key." has ".@$entries." columns with NULL values\n";
      foreach my $index(@$entries){
        print $index."\n";
      }
      print "\n\n";
      $count++;
      last if($count >= 10);
    }
  }
  my $has_trailing_hash = has_trailing_spaces($index_file);
  print "There are ".keys(%$has_trailing_hash)." lines in ".$index_file." with ".
      "trailing spaces on the end of columns\n";
  if($verbose){
    my $count;
    foreach my $key(keys(%$has_trailing_hash)){
      my $entries = $has_trailing_hash->{$key};
      print $key." has ".@$entries." columns with NULL values\n";
      foreach my $index(@$entries){
        print $index."\n";
      }
      print "\n\n";
    }
    $count++;
    last if($count >= 10);
  }
  
  my $column_check_results = check_column_syntax($index_file);
  print "There are ".keys(%$column_check_results)." lines with columns that have ".
      "syntax problems from ".$index_file."\n";
  if($verbose){
    my $count;
    my $name_to_index = index_to_name_hash();
    foreach my $key(keys(%$column_check_results)){
      print "Problem with some of the columns from ".$key."\n";
      my @indexes = keys(%{$column_check_results->{$key}});
      my @sorted = sort {$a <=> $b} @indexes;
    INDEX:foreach my $index(@sorted){
      my $name = $name_to_index->{$index};
      if(!$name){
	throw("Failed to get name from ".$index);
      }
      my $value = $column_check_results->{$key}->{$index};
      if(defined($value)){
        print $name." ".$index." has ".$value." rather that the correct form\n";
      }else{
        print $name." ".$index." column is undefined\n";
      }
    }
      $count++;
      last if($count >= 10);
    }
  }
  my $analysis_group_results = check_analysis_group($index_file);
  print "There are ".keys(%$analysis_group_results)." lines with problems for the ".
      "analysis group\n";
  if($verbose){
    my $count = 0;
    foreach my $key(keys(%$analysis_group_results)){
      print $key." ".$analysis_group_results->{$key}."\n";
      $count++;
      last if($count >= 10);
    }
  }
  my $population_results = check_population($index_file);
  print "There are ".keys(%$population_results)." lines with problems for the ".
      "population column\n";
  if($verbose){
    my $count = 0;
    foreach my $key(keys(%$population_results)){
      print $key." ".$population_results->{$key}."\n";
      $count++;
      last if($count >= 10);
    }
  } 
  #START LOG FILE WRITING
  #General syntax log: DATE TIME, Total number, Wrong column count, Undifined columns, Null strings columns, Trailing spaces on end, Syntax problems 
  open(DAT,">>log_syntax.csv") || die("Cannot Open File");
  print(DAT join(",",$theGMTime,"".keys(%$index_hash),"".keys(%$column_count_hash),"".keys(%$column_value_hash),"".keys(%$has_nulls_hash),"".keys(%$has_nulls_hash),"".keys(%$column_check_results)));
  close(DAT);
  
  if($log)
  {
    #start detailed log, like verbose only now into file	
    open(DAT,">log_syntax_ColumnCount_".$theGMTime.".csv") || die("Cannot Open File");
    foreach my $key(keys(%$column_count_hash)){ 
      print DAT $key." ".$column_count_hash->{$key}."\n";  
    }
    close(DAT);
    
    open(DAT,">log_syntax_UndifinedColumns_".$theGMTime.".csv") || die("Cannot Open File");	
    print DAT "Error\tKey\tIndex\n";
    if($column_value_hash && keys(%$column_value_hash)){
    KEY:foreach my $key(keys(%$column_value_hash)){
      my $values = $column_value_hash->{$key};
      foreach my $index(@$values){
        print DAT "Is undefined\t".$key."\t".$index."\n"
            if($index == 0 || $index == 1 || $index == 2 || $index == 3 ||
               $index == 5 || $index == 6 || $index == 10 || $index == 14 ||
               $index == 18);
      }
    }
    }
    close DAT;
    
    open(DAT,">log_syntax_NullStringsColumns_".$theGMTime.".csv") || die("Cannot Open File");
    print DAT "Error\tKey\tCSV Index\n";	
    foreach my $key(keys(%$has_nulls_hash)){
      my $entries = $has_nulls_hash->{$key};
      print DAT "Column with NULL values\t".$key."\t".join(",",@$entries);
      
    }
    close DAT;
    
    open(DAT,">log_syntax_HasTrailing_".$theGMTime.".csv") || die("Cannot Open File");
    print DAT "Error\tKey\tCSV Index\n";	
    foreach my $key(keys(%$has_trailing_hash)){
      my $entries = $has_trailing_hash->{$key};
      print DAT "Column with Trailings\t".$key."\t".join(",",@$entries);
    }	
    close DAT;
    
    open(DAT,">log_syntax_ColumnCheck_".$theGMTime.".csv") || die("Cannot Open File");
    print DAT "Error\tKey\tCSV Index\n";	
    my $name_to_index = index_to_name_hash();
    foreach my $key(keys(%$column_check_results)){
      
      
      my @indexes = keys(%{$column_check_results->{$key}});
      my @sorted = sort {$a <=> $b} @indexes;
    INDEX:foreach my $index(@sorted){
      my $name = $name_to_index->{$index};
      my $value = $column_check_results->{$key}->{$index};
      if(defined($value)){
        print DAT "Column with problem\t".$key."\t". $name." ".$index." has ".$value." rather that the correct form\n";
      			}else{
                          print DAT "Column with problem\t".$key."\t". $name." ".$index." column is undefined\n";
      			}
    }
    }
   }
}


sub useage{
  my ($args) = @_;
  print "Commandline was run_copy.pl ".join(" ", @$args)."\n";
  perldocs();
}

sub perldoc{
  exec('perldoc', $0);
  exit(0);
}



=pod

=head1 NAME

  OneKGenomes/scripts/meta/sequence_index_sanity.pl

=head1 SYNOPSIS

This script takes a sequence.index file as described in 
ftp://ftp.1000genomes.ebi.ac.k/vol1/ftp/README.sequence_data and does some sanity
checks on the contents

=head1 OPTIONS

    -index_file path to index file in above format
    -verbose binary flag for additional info about problems
    -syntax binary flag to indicate that the file syntax checks need to be run
    -filesystem binary flag to indicate that the filesystem checks need to be run
    -log binary flag for additional infor about problems written to log file
    -ftp_root root path for ftp site
    -ftp_path_must_contain string that the ftp paths must contain, by default it is
     sequence_read and the filesystem checks will only consider filepaths with
     sequence_read in them
    -dir_must_contain another level of regex to ensure only the correct files
     are compared. This was implemented to distiguish between ftp/data and 
     ftp/pilot_data/data
    -log writes detailed log files containing all errors and problems detected by -filesystem or/and -syntax
    -help print out the help

=head1 Sanity Checks

Filesystem

There is one filesystem check, compare_index_to_ftp. This generates 2 lists of paths
once based on the ftp_root and one based on the index file and compared them 
indicating which opaths are missing from the other list

Syntax

There are many syntax checks

check_column_count_sanity checks that the file contains the correct number of columns
has_empty_columns, check that the lines have defined values. This will count all 
undefined columns but when using the -verbose option will only compare about 
required values
has_null_strings, this will complain about columns which has the string NULL in them
has_trailing_spaces, this will complain about columns which have trailing spaces
check_column_syntax, this does a specific check for each column
0, fastq name, checks name is defined and  matches ^[S|E]RR.*\.fastq\.gz$
1, md5, checks is defined and is a 32 character long string
2, run_id, checks is defined and matches ^[S|E]RR\d+$
3, study_id, checks is defined and matches ^[S|E]RP\d+$
4, study name, no check
5, center name, checks is defined
6, submission_id checks is defined and matches ^[S|E]RA\d+$
7, submission_date checks if defined it matches ^\d{4}.\d{1,2}.\d{1,2}
8, sample_id checks if defined it matched ^[S|E]RS\d+$
9, sample_name, checks if defined and matches ^NA\d+$
10, population, no check
11, experiment_id, checks if defined and matches ^[S|E]RX\d+$
12, Instrument_platform, checks if defined
13, Instrument_model, no check
14, Library name, checks if defined
15, Run name, no check
16, Run Block Name, no check
17, Insert size, if defined check its just a number
18, Library Layout, check if defined and that is should be either SINGLE or PAIRED
19, Paired fastq file, if defined do check 0
20, Withdrawn, fail if not 1 or 0
21, Withdrawn date, no check
22, comment, no check
23, read count, must be a number or not available
24, base count, must be a number or not available
25, analysis group, must be one of these 2 strings low coverage, high coverage
or exon targetted

=head1 Examples

sequence_index_sanity.pl -index_file sequence.index -syntax -verbose
There are 17655 lines in /homes/laura/code/1000genomes/dcc/laura/march09_newsrr/docs/experiment_id_fix_sequence.index with undefined columns
There are 0 lines in /homes/laura/code/1000genomes/dcc/laura/march09_newsrr/docs/experiment_id_fix_sequence.index with NULL strings in columns
There are 0 lines in /homes/laura/code/1000genomes/dcc/laura/march09_newsrr/docs/experiment_id_fix_sequence.index with trailing spaces on the end of columns
There are 18 lines with columns that have syntax problems from /homes/laura/code/1000genomes/dcc/laura/march09_newsrr/docs/experiment_id_fix_sequence.index
Problem with some of the columns from data/NA19239/sequence_read/ERR001173_1.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form
Problem with some of the columns from data/NA19240/sequence_read/ERR000309_1.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form
Problem with some of the columns from data/NA19239/sequence_read/ERR001174_1.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form
Problem with some of the columns from data/NA19239/sequence_read/ERR001172_1.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form
Problem with some of the columns from data/NA19239/sequence_read/ERR001170_2.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form
Problem with some of the columns from data/NA19239/sequence_read/ERR001085_1.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form
Problem with some of the columns from data/NA19239/sequence_read/ERR001172_2.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form
Problem with some of the columns from data/NA19239/sequence_read/ERR000298_1.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form
Problem with some of the columns from data/NA19239/sequence_read/ERR001085_2.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form
Problem with some of the columns from data/NA19239/sequence_read/ERR001086_1.recal.fastq.gz
submission_date 7 has 2008-OCT-07 17:00:00 rather that the correct form

=cut

