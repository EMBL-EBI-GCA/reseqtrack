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

if ( !$input{collection} ) {
	throw("Please provide a collection name\n");
}

$input{host} = '1000genomes.ebi.ac.uk' if (!$input{host});
$input{update} = 0 if (!$input{update});
$input{collection_type} = "5PF_FASTQ" if (!$input{collection_type});
$input{output_file_type} = "SAM_CHUNK" if (!$input{output_file_type});

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd`;
	chomp $input{output_dir};
}

my @tmp = split(/\_/, $input{collection});
my $pop = $tmp[0];
my $sample_cnt = $tmp[1];
my $outdir_pop = $input{output_dir} . "/$pop";
mkpath($outdir_pop) unless (-e $outdir_pop);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
);

my $fastq_col = $db->get_CollectionAdaptor->fetch_by_name_and_type($input{collection}, $input{collection_type});
my $fastq_files = $fastq_col->others;

$db->dbc->disconnect_when_inactive(1); 

foreach my $fastq ( @$fastq_files ) {
	
	throw("Chunky fastq file " . $fastq->name . " does not exist") unless (-e $fastq->name);
	my $sam_out = $fastq->name . ".sam";
	throw("chunky sam file $sam_out already exist") if (-e $sam_out);
	
	#my $stampy_bin = "/nfs/1000g-work/G1K/work/bin/stampy/stampy.py";
	my $stampy_bin = "/nfs/production/reseq-info/work/bin/stampy/stampy.py";
	my $stampy_hash_stub = "/nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37";
	my $command = "$stampy_bin -g $stampy_hash_stub -h $stampy_hash_stub --norefoutput --inputformat=fasta -M " . $fastq->name . " -o $sam_out";
	print "Run stampy command if -run:\n$command\n";
	`$command` if $input{run};
	my $exit = $?>>8;
	throw("Run stampy failed for $fastq\n") if ($exit >=1);	
	
	if (-e $sam_out) {
		my $sam_objs = save_file($sam_out, $input{output_file_type});
		my @tmp3 = split(/\./, $input{collection});
		my $total_file_cnt = pop @tmp3;
		my $collection_name = $tmp3[0] . "_$total_file_cnt" . ".sam";
		print "SAM_CHUNK collection name is $collection_name\n";
		my $collection =  ReseqTrack::Collection->new(
			-name => $collection_name,
			-others => $sam_objs,
			-type => $input{output_file_type},  
			-table_name => 'file',
			);
		save_collection($collection, $collection_name, $input{output_file_type});		
	}
	else {
		throw("chunky sam file $sam_out does not exist");
	}	
}	

############ SUBS ##########
sub save_collection { ### FIXME, need more work!
	my ($co, $c_name, $c_type) = @_;
	#print "collection to store is $c_name, type $c_type\n";
	my $exist = $db->get_CollectionAdaptor->fetch_by_name_and_type($c_name, $c_type);
	if ( !$exist ) {
		$db->get_CollectionAdaptor->store($co) if ($input{save_collection});  ## new files will be added in without over write the old ones.
		print "collection $c_name, type $c_type is stored\n";
	}
	else{
		$db->get_CollectionAdaptor->store($co, 1) if ($input{save_collection});
		print "collection $c_name, type $c_type is updated\n";
	}	
	return 1;
}			


sub save_file {
	my ($file, $type) = @_;
	my $f_objs;		
	my $loader = ReseqTrack::Tools::Loader::File->new
		  (	-file => [$file],
			-do_md5 => 1,
			-hostname => $input{host},
			-db => $db,
			-assign_types => 0,
			-check_types => 0,
			-type => $type,
			-update_existing => $input{update},
		) if ( $input{store} );
			
	if($input{store}){
		$loader->process_input();
		$loader->create_objects();
		$loader->sanity_check_objects();
		$f_objs = $loader->load_objects();
	}
	return $f_objs;
}
	
=pod
perl $ZHENG_RB_VC/scripts/process/cortex_run_stampy.pl $WRITE_DB_ARGS -dbname zheng_run_cortex_test\
-output_dir /nfs/1000g-work/G1K/work/zheng/cortex/post_process \
-collection_type 5PF_FASTQ \
-output_file_type SAM_CHUNK \
-collection PEL_50_bubbled_sample_ctxs.1.genotype.5pflank
