#!/sw/arch/bin/perl -w

use strict;

use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::RunMetaInfoUtils;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Tools::FileUtils;
use Getopt::Long;
use File::Basename;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $type;
my $table_name = "file";
my $output_file;
my $root_trim = '/nfs/1000g-archive/vol1/ftp/';
my $single_run_id;
my @skip_study_ids;
my $study_collection_type = 'STUDY_TYPE';
my @print_status;
my $help;

&GetOptions(
	    'dbhost=s'     		=> \$dbhost,
	    'dbname=s'      	=> \$dbname,
	    'dbuser=s'      	=> \$dbuser,
	    'dbpass=s'      	=> \$dbpass,
	    'dbport=s'      	=> \$dbport,
	    'type=s' 			=> \$type,
	    'table_name=s' 		=> \$table_name,
	    'output_file=s' 	=> \$output_file,
	    'root_trim=s' 		=> \$root_trim,
	    'skip_study_id=s@' 	=> \@skip_study_ids,
	    'run_id=s' 			=> \$single_run_id,
	    'study_collection_type:s' => \$study_collection_type,
        'print_status=s'  	=> \@print_status,
       	'help!' 			=> \$help,
);

if($help){
  useage();
}

if (!@print_status) {
  push(@print_status, 'public');
}
  
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );
    
my %skip_study_id;
foreach my $study_id(@skip_study_ids){
  $skip_study_id{$study_id} = 1;
}

my $ca = $db->get_CollectionAdaptor;
my $study_collections = $ca->fetch_by_type($study_collection_type);
my %study_collection_hash;
foreach my $collection(@$study_collections){
	next if ($collection->name 	!~ /high coverage cgi/ ); # skip runs that are not Complete Genomics runs
	my $others = $collection->others;
	foreach my $other(@$others){
		$study_collection_hash{$other->run_id} = $collection->name;
	}
}

my $rmi_a = $db->get_RunMetaInfoAdaptor;
my $meta_infos;
if($single_run_id){
  my $rmi = $rmi_a->fetch_by_run_id($single_run_id);
  $meta_infos = [$rmi];
}else{
  $meta_infos = $rmi_a->fetch_all;
}
my @sorted = sort{$a->run_id cmp $b->run_id} @$meta_infos;

my $fa = $db->get_FileAdaptor;
my %file_hash;
if($table_name eq 'file'){
  my $files = $fa->fetch_by_type($type);
  foreach my $file(@$files){
    my $dir = dirname( $file->name );
    my ($sample) = split (/\_/, $file->filename);
    #print "sample is $sample, dir is $dir\n";
    $file_hash{$sample} = $dir;
  }
}

open (OUT, ">", $output_file) || throw("Cannot open output file $output_file");

my $methods = index_method_array();
my $headers = index_method_array();
my @header;
foreach my $column ( @$headers ) {
	$column =~ tr/a-z/A-Z/;
	push @header, $column;
}	
print OUT join ("\t", @header) . "\t";
print OUT "DIRECTORY\tANALYSIS_GROUP\n";
  
my %index_lines;
foreach my $meta_info(@sorted){

	next if ( !$study_collection_hash{$meta_info->run_id} );

	if($single_run_id){
     #print STDERR "Comparing ".$meta_info->run_id." to ".$single_run_id."\n";
		next unless($meta_info->run_id eq $single_run_id);
	}
	
	if(keys(%skip_study_id)){
		next if($skip_study_id{$meta_info->study_id});
	}

	next unless (grep {$_ eq $meta_info->status} @print_status); ## This is to skip the runs with status of non "public"
	
	next if ($meta_info->center_name !~ /CompleteGenomic/i);

	my $string;

	foreach my $method(@$methods){ # these are column headers (or column names for the run_meta_info table)  
		$string .= $meta_info->$method."\t";
	}
  
	my $cg_dir_path;
	if ( $file_hash{$meta_info->sample_name} ) {
		$cg_dir_path = $file_hash{$meta_info->sample_name};
	}
	elsif ( $file_hash{$meta_info->sample_name . "-1"} ) {
		$cg_dir_path = $file_hash{$meta_info->sample_name . "-1"};
	}  
	elsif ( $file_hash{$meta_info->sample_name . "-C1"} ) {
		$cg_dir_path = $file_hash{$meta_info->sample_name . "-C1"};
	}
	else {
		throw("No data directory is found for run" . $meta_info->run_id);
	}   
	$cg_dir_path =~ s/$root_trim//;
	$string .= $cg_dir_path . "\t"; 
	$string .= $study_collection_hash{$meta_info->run_id} . "\t";
	print OUT "$string\n";
}

=pod
=head1 NAME

ReseqTrack/scripts/run_meta_info/dump_CompleteGenomics_sequence_index.pl

=head1 SYNOPSIS

This script dumps a CompleteGenomics.sequence.index file based on the run meta information information 
in the given database conntecting to Complete Genomics data directory based on user-input file type

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

-root_trim, this is the root path which needs to be trimmed

-single_run_id, this is to allow a single runs index to be dumped which is useful
for debugging purposes

-print_status, lines will be printed if the run_meta_info status is equal to this.
  Default is 'public'.  Can be specified multiple times for multiple acceptable statuses.

-help, binary flag to get the perldocs printed

=head1 Examples

$ perl $ZHENG_RT/scripts/run_meta_info/dump_CompleteGenomics_sequence_index.pl $WRITE_DB_ARGS -dbname g1k_archive_staging_track -table_name file -type CG_BAM -output_file 20130129.CompleteGenomics.sequence.index

