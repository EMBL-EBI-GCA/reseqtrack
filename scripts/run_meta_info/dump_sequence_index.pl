#!/sw/arch/bin/perl -w

use strict;

use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::RunMetaInfoUtils;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Tools::FileUtils;
use Getopt::Long;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $type;
my $table_name;
my $help;
my $output_file;
my $trim_paths = 1;
my $root_trim = '/nfs/1000g-archive/vol1/ftp/';
my $single_run_id;
my @skip_study_ids;
my $test_id;

&GetOptions( 
  'dbhost=s'      => \$dbhost,
  'dbname=s'      => \$dbname,
  'dbuser=s'      => \$dbuser,
  'dbpass=s'      => \$dbpass,
  'dbport=s'      => \$dbport,
  'type=s' => \$type,
  'table_name=s' => \$table_name,
  'help!' => \$help,
  'output_file=s' => \$output_file,
  'trim_paths!' => \$trim_paths,
  'root_trim=s' => \$root_trim,
  'skip_study_id=s@' => \@skip_study_ids,
  'run_id=s' => \$single_run_id,
    );

if($help){
  useage();
}
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );


my %skip_study_id;
my %staging_area;
foreach my $study_id(@skip_study_ids){
  $skip_study_id{$study_id} = 1;
}
#get meta info
#find files or collections based on type
#create index lines if there are or aren't files
#print lines
#throw("Don't have a single run id") unless($single_run_id);
my $rmi_a = $db->get_RunMetaInfoAdaptor;
my $meta_infos = $rmi_a->fetch_all;
my @sorted = sort{$a->run_id cmp $b->run_id} @$meta_infos;

my $fa = $db->get_FileAdaptor;
my %file_hash;
if($table_name eq 'file'){
  my $files = $fa->fetch_by_type($type);
  foreach my $file(@$files){
    my ($run_id) = $file->filename =~ /([E|S]RR\d+)/;
    push(@{$file_hash{$run_id}}, $file);
  }
}
my %index_lines;
 META_INFO:foreach my $meta_info(@sorted){
   print "Have ".$meta_info->run_id."\n";
   my $analysis_group = get_analysis_group($meta_info);
   if($single_run_id){
     #print STDERR "Comparing ".$meta_info->run_id." to ".$single_run_id."\n";
     next META_INFO unless($meta_info->run_id eq $single_run_id);
   }
   print "Running dump on ".$meta_info->run_id."\n";
   if(keys(%skip_study_id)){
     next META_INFO if($skip_study_id{$meta_info->study_id});
   }
   $index_lines{$meta_info->run_id} = [] unless($index_lines{$meta_info->run_id});
   if($meta_info->status eq 'suppressed' || $meta_info->status eq 'cancelled'){
      my $line = create_suppressed_index_line($meta_info, undef, $analysis_group);
     #print $line."\n";
     push(@{$index_lines{$meta_info->run_id}}, $line);
   }elsif($meta_info->status eq 'public'){
     my $files;
      if($table_name eq 'file'){
       $files = $file_hash{$meta_info->run_id};
       my @tmp;
       foreach my $file(@$files){
	 if($file->name =~ /$root_trim/){
	   push(@tmp, $file);
	 }else{
	   $staging_area{$meta_info->run_id} = 1;
	 }
       }
       $files = \@tmp;
     }elsif($table_name eq 'collection'){
       throw("Why are you using collections, this is a bad idea");
       my $collections = get_file_collections_associated_with_run($meta_info, $type);
       foreach my $collection(@$collections){
         push(@$files, @{$collection->others});
       }
     }
     
     #### This bit is to determine what comment to put in sequence index file when a line is obseleted
     my $eca = $db -> get_EventCompleteAdaptor;
     my $all_collections = get_file_collections_associated_with_run($meta_info); #get all collection regardless of type
     my $complete; 
   
     foreach my $coll (@$all_collections) {
       #$complete = $eca -> fetch_by_object_name_and_event($coll->name, $event);
       if ($coll->type eq "ARCHIVE_FASTQ") {
         $complete = $eca -> fetch_by_other_id($coll->dbID, "collection");
         #print "run id: " . $meta_info->run_id . " Complete event for collection " . $coll->dbID . " " . $coll->type . " is $complete deref:" . @$complete . "\n";
       }
     }
     

     if(!$files || @$files == 0){
       my $warning = $meta_info->run_id." seems to have no files associated with it";
       $warning .= " for type ".$type if($type);
       warning($warning);
      my $line;
       if ( !$all_collections || @$all_collections ==0 ) {
         $line = create_suppressed_index_line($meta_info, 'NOT YET AVAILABLE FROM ARCHIVE', $analysis_group);
       }elsif(@$all_collections == 1 && (!$complete || @$complete == 0)){
         $line = create_suppressed_index_line($meta_info, 'NOT YET AVAILABLE FROM ARCHIVE', $analysis_group);
       }elsif($staging_area{$meta_info->run_id}){
	 $line = create_suppressed_index_line($meta_info, 'NOT YET AVAILABLE FROM ARCHIVE', $analysis_group);
       }elsif(@$all_collections >= 2 && (!$complete || @$complete == 0)){
         $line = create_suppressed_index_line($meta_info, 'FAILED ONE OF THE DCC PROCESSES', $analysis_group);
       }elsif(@$all_collections >= 2 && @$complete == 1){
         $line = create_suppressed_index_line($meta_info, 'FAILED GENOTYPE QC', $analysis_group);
       }elsif(@$all_collections == 1 && @$complete ==1){
         $line = create_suppressed_index_line($meta_info, 'FAILED ONE OF THE DCC PROCESSES', $analysis_group);
       }else{
         throw("Have ".@$all_collections." collections and ".@$complete." ".
               "event complete objects associated with ".$meta_info->run_id
               ." Not sure what sort of line to print")
       }
       push(@{$index_lines{$meta_info->run_id}}, $line);
     }else{
       #print STDERR "Have ".@$files." for ".$meta_info->run_id."\n";
       my ($mate1, $mate2, $frag) = assign_file_objects($files);
       #print STDERR "Mate 1 ".$mate1->name."\n" if($mate1);
       #print STDERR "Mate 2 ".$mate2->name."\n" if($mate2);
       #print STDERR "Frag ".$frag->name."\n" if($frag);
       if($frag){
         my ($read_count, $base_count) = get_count_stats($frag);
         my $line = create_index_line($frag->name, $frag->md5, $meta_info, undef, 0,
                                      undef, undef, $read_count, $base_count, 
                                      $analysis_group);
         push(@{$index_lines{$meta_info->run_id}}, $line);
       }
       if($meta_info->library_layout eq 'SINGLE'){
        if($mate1 || $mate2){
          warning("Can't handle a single ended library ".$meta_info->run_id.
                  " with either ".$mate1." ".$mate2." defined");
          next META_INFO;
        }
       }elsif($meta_info->library_layout eq 'PAIRED'){
        if($mate1 && $mate2){
          my ($mate1_read_count, $mate1_base_count) = get_count_stats($mate1);
          my $mate1_line = create_index_line($mate1->name, $mate1->md5, $meta_info, 
                                             $mate2->name, 0, undef, undef, 
                                             $mate1_read_count, $mate1_base_count,
                                             $analysis_group);
          my ($mate2_read_count, $mate2_base_count) = get_count_stats($mate2);
          my $mate2_line = create_index_line($mate2->name, $mate2->md5, $meta_info, 
                                             $mate1->name, 0, undef, undef,
                                             $mate2_read_count, $mate2_base_count,
                                             $analysis_group);
          push(@{$index_lines{$meta_info->run_id}}, ($mate1_line, $mate2_line));
        }else{
          warning("Dont have mate1 and/or mate2 files for ".$meta_info->run_id);
        }
      }else{
        throw("Don't know what to expect for ".$meta_info->library_layout);
      }
     }
   }
}

my $fh = \*STDOUT;

if($output_file){
  open(FH, ">".$output_file) or throw("Failed to open ".$output_file." $!");
  $fh = \*FH;
}

my $header = return_header_string();
print $fh $header;

foreach my $meta_info(@sorted){
  my $lines = $index_lines{$meta_info->run_id};
  if($lines && @$lines > 0){
    foreach my $line(@$lines){
      if($trim_paths){
        my @values = split /\t/, $line;
        $values[0] =~ s/$root_trim//;
        $values[19]  =~ s/$root_trim// if($values[19]);
        print $fh join("\t", @values)."\n";
      }else{
        print $fh $line."\n";
      }
    }
  }else{
    print STDERR "Have no lines for ".$meta_info->run_id."\n" unless($single_run_id);
  }
}


close($fh);

=pod

=head1 NAME

ReseqTrack/scripts/run_meta_info/dump_sequence_index.pl

=head1 SYNOPSIS

This script dumps a sequence.index file based on the run meta information information 
in the given database conntecting to files based on types from either the file or the 
collection table

=head1 OPTIONS

-dbhost, the name of the mysql-host

-dbname, the name of the mysql database

-dbuser, the name of the mysql user

-dbpass, the database password if appropriate

-dbport, the port the mysql instance is running on, this defaults to 4197 the 
    standard port for mysql-g1kdcc.ebi.ac.uk

-table_name, this is to indicate if you want to associate your runs with files or collections

-type, this is the type string you will use to fetch either files or collections

-skip_study_id, this is a option to allow specific studies to be skipped, it can
appear multiple times on the commandline

-output_file, this is where the sequence.index file is dumped out to

-trim_path, this indicates that a root path should be trimmed from the file paths

-root_trim, this is the root path which needs to be trimmed

-single_run_id, this is to allow a single runs index to be dumped which is useful
for debugging purposes

-help, binary flag to get the perldocs printed

=head1 Examples


perl ReseqTrack/scripts/run_meta_info/dump_sequence_index.pl -dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database -table_name file -type FILTERED_FASTQ


=head1 Other useful scripts

ReseqTrack/scripts/run_meta_info/dump_sequence_index_stats.pl

=cut
