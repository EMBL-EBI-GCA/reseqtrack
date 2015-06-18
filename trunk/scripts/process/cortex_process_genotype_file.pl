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
use ReseqTrack::Tools::RunSplit

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
	
	'action=s',
	'collection=s',
	'collection_type=s',
	
	'line_cnt_per_fastq=i',
	'split_fastq_program=s',
	
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
$input{collection_type} = "GENOTYPE" if (!$input{collection_type});
	
if ( !$input{collection}  ) {
	throw("Please provide a collection name\n");
}
if ( !$input{action}  ) {
	throw("Please provide an action, either create_aux_cov_file_fr_genotype or make_and_split_5pf_fastq_fr_genotype\n");
}

if ($input{action} eq "create_aux_cov_file_fr_genotype") {
	$input{output_file_type} = "COVG_FILE" if ( !$input{output_file_type} );
}
elsif (	$input{action} eq "make_and_split_5pf_fastq_fr_genotype") {
	#$input{split_fastq_program} = "/nfs/1000g-work/G1K/work/zheng/reseqtrack/c_code/split/split" if (!$input{split_fastq_program});	
	$input{split_fastq_program} = "/nfs/production/reseq-info/work/zheng/reseqtrack/trunk/c_code/split/split" if (!$input{split_fastq_program});
	$input{line_cnt_per_fastq} = 1000000 if ( !$input{line_cnt_per_fastq}); ## split 5pf fastq file intp 1M-line chunks
	$input{output_file_type} = "5PF_FASTQ" if ( !$input{output_file_type} );
}

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd`;
	chomp $input{output_dir};
}

my @tmp = split(/\_/, $input{collection});
my $outdir_pop = $input{output_dir} . "/$tmp[0]";
mkpath($outdir_pop) unless (-e $outdir_pop);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
);

my $col_obj = $db->get_CollectionAdaptor->fetch_by_name_and_type($input{collection}, $input{collection_type});
my $file_objs = $col_obj->others;

$db->dbc->disconnect_when_inactive(1); 

if ( $input{action} eq "create_aux_cov_file_fr_genotype" ) {
	foreach my $file_obj ( @$file_objs ) {
		throw("Genotype file " . $file_obj->name . " does not exist") unless (-e $file_obj->name);
		print "genotype file is " . $file_obj->name ."\n";
		my $covg_file = create_aux_covg_file($file_obj);
		save_file($covg_file, "COVG_FILE") if ($input{run});		
	}	
}
elsif (  $input{action} eq "make_and_split_5pf_fastq_fr_genotype"  ) {
	foreach my $file_obj ( @$file_objs ) {
		
		throw("Genotype file " . $file_obj->name . " does not exist") unless (-e $file_obj->name);
		print "genotype file is " . $file_obj->name ."\n";
	
		my $flank_seq = $outdir_pop . "/" . $input{collection} . ".5pf.fastq"; 
		#my $command = "perl /nfs/1000g-work/G1K/work/bin/cortex/scripts/analyse_variants/make_5p_flank_file.pl " . $file_obj->name . " > $flank_seq";
		my $command = "perl /nfs/production/reseq-info/work/bin/cortex/scripts/analyse_variants/make_5p_flank_file.pl " . $file_obj->name . " > $flank_seq";
		print "Making 5pflanking seq fastq file..... if -run\n";
		print "$command\n";
		`$command` if $input{run};
		my $exit = $?>>8;
		throw("Making 5p flanking fastq failed for " . $file_obj->name . "\n") if ($exit >=1);	
			
		if (-e $flank_seq ) {
			my $split_obj = split_fastq($flank_seq, $outdir_pop);
			my $output_file_hash = $split_obj->output_file_hash;
			foreach my $ori_fastq ( keys %$output_file_hash ) {
				my $split_file_cnt = keys ( %{$output_file_hash->{$ori_fastq}} );
				foreach my $i (keys %{$output_file_hash->{$ori_fastq}} )  {
					print "split file is $output_file_hash->{$ori_fastq}->{$i}\n";
					my $fobjs = save_file($output_file_hash->{$ori_fastq}->{$i}, "5PF_FASTQ") if ($input{run}); ### stampy is ok with fastq.gz, in addition to fastq
					my $collection_name = basename($output_file_hash->{$ori_fastq}->{$i}) . ".total$split_file_cnt";
					print "5PF_FASTQ collection name is $collection_name\n";
					my $collection =  ReseqTrack::Collection->new(
						-name => $collection_name,
						-others => $fobjs,
						-type => "5PF_FASTQ",  
						-table_name => 'file',
					);
					save_collection($collection, $collection_name,"5PF_FASTQ") if ($input{run});	
				}		
			}
			`rm $flank_seq` if ($input{run});	
		}
		else {
			throw("Unsplit 5PF fastq file $flank_seq does not exist\n");
		}
	}			
}
	
##### SUBS ####
sub create_aux_covg_file {
	my ($gf) = @_;
	my $covg_for_classifier_file = $gf->name . ".covg_for_classifier"; ## That's where Zam's code will write the output file
	throw("coverage file $covg_for_classifier_file already exist, delete it before run") if (-e $covg_for_classifier_file);
	my $gf_basename = basename($gf->name);

	my @tmp = split(/\_/, $gf_basename);
	my $sample_number_plus_1 = $tmp[1] + 1; 

	#my $command2 = "perl /nfs/1000g-work/G1K/work/bin/cortex/scripts/analyse_variants/make_covg_file.pl " . $gf->name . " $sample_number_plus_1 0";
	my $command2 = "perl /nfs/production/reseq-info/work/bin/cortex/scripts/analyse_variants/make_covg_file.pl " . $gf->name . " $sample_number_plus_1 0";

	eval {
		print "Running \n$command2 if -run\n";
		`$command2` if ($input{run});	
	};
	throw("make_covg_file.pl failed for " . $gf->name . ", $@\n") if $@;
		
	return $covg_for_classifier_file;
}

sub split_fastq {
	my ($fastq, $od) = @_;
	my %line_count_hash;
	$line_count_hash{$fastq} = $input{line_cnt_per_fastq}; # split the fastq file into sub files, 1M lines each if the input is 1M
	my $run_split = ReseqTrack::Tools::RunSplit->new(
		-program => $input{split_fastq_program},
		-working_dir => $od,
		-line_count_hash => \%line_count_hash,
		-job_name => "NA",
		-working_dir => $od );
	print "Splitting the 5pflank fastq into 1M lines chunks if -run.......\n";
	$run_split->run if ($input{run});				
	return $run_split;		
}

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
perl $ZHENG_RB_VC/scripts/process/cortex_process_genotype_file.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-output_dir /nfs/1000g-work/G1K/work/zheng/cortex/post_process \
-collection_type GENOTYPE \
-action create_aux_cov_file_fr_genotype \
-output_file_type COVG_FILE \
-collection ACB_64_bubbled_sample_ctxs.genotype \
-run \
-store \
-save_collection \
-update

perl $ZHENG_RB_VC/scripts/process/cortex_process_genotype_file.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-output_dir /nfs/1000g-work/G1K/work/zheng/cortex/post_process \
-collection_type GENOTYPE \
-action make_and_split_5pf_fastq_fr_genotype \
-output_file_type 5PF_FASTQ \
-collection ACB_64_bubbled_sample_ctxs.genotype \
-run \
-store \
-save_collection \
-update
	
