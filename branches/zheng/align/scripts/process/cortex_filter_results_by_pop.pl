#!/usr/bin/env perl

use strict;
use warnings;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::GeneralUtils;
use Getopt::Long;
use File::Basename;
use File::Path;
use ReseqTrack::Tools::Loader::File;
use Time::Local;
use ReseqTrack::Tools::RunVariantCall;


my $start_time = time();

$| = 1;

my (
	%input,
);		


GetOptions(
	\%input,    
	'dbhost=s',   
	'dbname=s',      
	'dbuser=s',
	'dbpass=s', 
	'dbport=s',

	'collection=s',
	'collection_type=s',
	
  	'output_dir=s',
  	'output_file_type:s',

  	'store!',
  	'save_collection!',
 	'update!',
 	'host:s',
  	'run!',		
);

$input{host} = '1000genomes.ebi.ac.uk' if (!$input{host});
$input{update} = 0 if (!$input{update});
$input{collection_type} = "AUX" if (!$input{collection_type});
$input{output_file_type} = "CLASSIFIED_CHUNK" if (!$input{output_file_type});

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd`;
	chomp $input{output_dir};
}
if ( !$input{collection} ) {
	throw("Please provide a collection name\n");
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
);

my $aux_col = $db->get_CollectionAdaptor->fetch_by_name_and_type($input{collection}, $input{collection_type});
my $aux_files = $aux_col->others;

my $table_file;
my $covg_file;
foreach my $aux_file ( @$aux_files ) {
	$table_file = $aux_file if ($aux_file->type eq "TABLE_FILE");
	$covg_file = $aux_file if ($aux_file->type eq "COVG_FILE");
}

my @tmp = split(/\_/, $aux_col->name);
my $pop = $tmp[0];
my $total_var_cnt = $tmp[1];
my $colour_cnt = $tmp[2];
my $start = $tmp[3];

$db->dbc->disconnect_when_inactive(1); 

my $odir = $input{output_dir} . "/$pop";
mkpath($odir) unless (-e $odir);
my $output_file = $odir . "/" . $aux_col->name . "_genotype.classified";

throw("Chunky classified var file $output_file already exist") if (-e $output_file);

my $command = 	"cat /nfs/1000g-work/G1K/work/bin/cortex/scripts/analyse_variants/classifier.parallel.ploidy_aware.R | ";
$command .= "R --vanilla --args $start 500000 " . $covg_file->name . " $total_var_cnt $colour_cnt 1 ";
#$command .= "R --vanilla --args $start 1000 " . $covg_file->name . " $total_var_cnt $colour_cnt 1 ";
$command .= $table_file->name . " 3000000000 31 2 ";
$command .= "$output_file";

print $command . "\n";

`$command` if ($input{run});
my $exit = $?>>8;
throw("filter variant calls by population failed for " . $aux_col->name . "\n") if ($exit >=1);	

if ( -e $output_file ) {
	my $loader = ReseqTrack::Tools::Loader::File->new
			(
			   -file => [$output_file],
			   -do_md5 => 1,
			   -size => (-s $output_file),
			   -hostname => $input{host},
			   -db => $db,
			   -assign_types => 0,
			   -check_types => 0,
			   -type => $input{output_file_type},
			   -update_existing => $input{update},
	);
			
	my $output_file_objects;
	if($input{store}){
		  $loader->process_input();
		  $loader->create_objects();
		  $loader->sanity_check_objects();
		  $output_file_objects = $loader->load_objects();
	}	
		
	if ( $input{save_collection} ) {
		my $collection_name = $pop . "_$total_var_cnt". "_$colour_cnt" . ".classified";
		my $collection =  ReseqTrack::Collection->new(
			-name => $collection_name,
			-others => $output_file_objects,
			-type => $input{output_file_type},  
			-table_name => 'file',
		);
		$db->get_CollectionAdaptor->store($collection, 1);  ## new files will be added in without over write the old ones.
	}
}	
	
=pod
perl $ZHENG_RB_VC/scripts/process/cortex_filter_results_by_pop.pl $WRITE_DB_ARGS -dbname zheng_run_cortex -collection_type AUX -output_file_type CLASSIFIED_CHUNK -collection ACB_1688175_65_1 -output_dir /tmp -run -store 
-save_collection -update




