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

use ReseqTrack::Tools::RunVariantCall::CallByCortex;

my $start_time = time();

$| = 1;

my (
	%input,
	$se_list,
	$pe1_list,
	$pe2_list,
	$pop,
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
	
	'ref_binary:s',
	
	'bubble_fasta:s',
	
	'colour_list:s', ## name has to contain the word "colour_list"
	
	'executable=s',
	'kmer_size=i',
	'mem_height=i',
	'mem_width=i',
	'dump_binary!',  
	'dump_covg_distribution!', 
	'detect_bubbles1:s', 
	'parameters=s%',
	
	'seq_index=s',
	
  	'output_dir=s',
  	'output_file_type:s',
  	'tabix_dir=s',

  	'store!',
  	'save_collection!',
 	'update!',
 	'host:s',
  	'help!',		
);

$input{update} = 0 if (!$input{update});

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd`;
	chomp $input{output_dir};
}
if ( $input{collection} && ! $input{seq_index}) {
	throw("Please provide sequence index, so the output can be organized by populations\n");
}
	
$input{host} = '1000genomes.ebi.ac.uk' if (!$input{host});

if (!$input{output_file_type} && $input{store} ) {
	throw("Please provide a file type if you want to store the output file in the database");
}	

if (!$input{tabix_dir}) {
	#$input{tabix_dir}= "/nfs/1000g-work/G1K/work/bin/tabix/";
	$input{tabix_dir}= "/nfs/production/reseq-info/work/bin/tabix/";
}	

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
);

my @inputs;
my ($samples_to_look_for, $sample_to_pop_mapping) = parse_seq_index($input{seq_index});
if ( 	$input{collection} && $input{collection_type} ) {
	my $collection =
	    $db->get_CollectionAdaptor->fetch_by_name_and_type($input{collection},
	                                                        $input{collection_type});
	throw("No collection found for $input{collection} and $input{collection_type}") if (!$collection); 
	
	foreach my $list (@{$collection->others}){
	    print "list: " . $list->name . "\n";
	    push @inputs, $list->name;
	}
	my @tmp = split(/\./, $input{collection});
	my $possible_sam_name = $tmp[0]; 
	if ($sample_to_pop_mapping->{$possible_sam_name} ) {
		$pop = $sample_to_pop_mapping->{$possible_sam_name};
	}
	else {
		my @tmp2 = split(/\_/, $input{collection});
		$pop = $tmp2[0] if ($input{collection} !~ /ref/);
	}	
	#$pop = $sample_to_pop_mapping->{$input{collection}}; #### FIXME:  BUG, this won't get pop as most collection name is not a sample name
}
elsif ( $input{colour_list} ) {
	push @inputs, $input{colour_list};
	if ($input{colour_list} =~ /^NA|^HG/) {
		my @tmp = split (/\./, $input{colour_list});
		my $sample = $tmp[0]; 
		$pop = $sample_to_pop_mapping->{$sample};
	}	
}	
elsif ( $input{bubble_fasta} ) {
	throw("Please name the bubble fasta file xxxx.se_list") if ($input{bubble_fasta} !~ /se_list/i);
	push @inputs, $input{bubble_fasta};
}	

$db->dbc->disconnect_when_inactive(1); 

### FIXME, define more ways to take in different input files

my $variant_calling_obj = ReseqTrack::Tools::RunVariantCall::CallByCortex->new(
	-program					=> $input{program},
	-executable					=> $input{executable},
	-input_files				=> \@inputs,
	-kmer_size 					=> $input{kmer_size},	
	-mem_height					=> $input{mem_height},
	-mem_width					=> $input{mem_width},
	-parameters					=> $input{parameters},
	-dump_binary				=> $input{dump_binary},
	-dump_covg_distribution		=> $input{dump_covg_distribution},
	-detect_bubbles1			=> $input{detect_bubbles1},
	-sample						=> $input{collection},
	-working_dir				=> $input{output_dir},
	-save_files_from_deletion	=> $input{save_files_from_deletion},
	-collection					=> $input{collection},
	-collection_type			=> $input{collection_type},
	-ref_binary					=> $input{ref_binary},
	-population					=> $pop,
);

$variant_calling_obj->run;

my $outfiles = $variant_calling_obj->output_files;							
							
#$db->dbc->disconnect_when_inactive(0);

foreach my $outfile ( @$outfiles ) {  
	
	if ( -e $outfile ) {	
		
		if ($outfile =~ /vcf$|vcf.gz$/ ) {
			my $zipped_vcf = bgzip_and_index($outfile, $input{tabix_dir});
			$outfile = $zipped_vcf;
		}		
			
		my $loader = ReseqTrack::Tools::Loader::File->new
		  (
		   -file => [$outfile],
		   -do_md5 => 1,
		   -hostname => $input{host},
		   -db => $db,
		   -assign_types => 0,
		   -check_types => 0,
		   -type => $input{output_file_type},
		   -update_existing => $input{update},
		  ) if ( $input{store} );
		
		my $file_objects;
		if($input{store}){
		  $loader->process_input();
		  $loader->create_objects();
		  $loader->sanity_check_objects();
		  $file_objects = $loader->load_objects();
		}

		my $collection_name;
		if ($input{save_collection} ) {  
			$collection_name = basename($outfile);
			my $collection =  ReseqTrack::Collection->new(
				  -name => $collection_name,
				  -others => $file_objects,
				  -type => $input{output_file_type},  
				  -table_name => 'file',
				);
			print "collection to store is $collection_name\n";
			$db->get_CollectionAdaptor->store($collection);  ## new files will be added in without over write the old ones.
			#### FIXME, dictate which process should store collection, which should not.
		}	
		
		print "*********************************************************************************************************************************\n";
		if ($input{store}) {
			print "**** Output file $outfile has been stored in the database\n";
		}
		else {
			print "**** Output file path is $outfile; it is NOT stored in the database as -store is not specified\n";
		}	
		if ( $input{save_collection}  ) {
		    print "**** Output file " . basename($outfile) . " has been stored in collection table as $collection_name\n";
		 }     
		print "*********************************************************************************************************************************\n";  
	}	
	else {
		throw("Output file " . $outfile . " does not exist");
	}	
}

my $end_time = time();
print "Job took " . ($end_time - $start_time)/60 . " minutes\n";
					
#### SUBS #####

sub bgzip_and_index {
	my ($vcf, $tabix_program) = @_;
	
	$tabix_program =~ s/\/$//;
	my $tabix = $tabix_program . "/tabix";
	my $bgzip = $tabix_program . "/bgzip";
	
	my $zip_vcf;
	if ($vcf !~ /gz$/) {
		`$bgzip -f $vcf`;
		my $exit = $?>>8;
		throw("bgzip failed for $vcf\n") if ($exit >=1);
		$zip_vcf = $vcf . ".gz";
	}
	else {
		$zip_vcf = $vcf;
	}
	
	my $index_file = $vcf . ".tbi";
	eval {
		`$tabix -p vcf $zip_vcf` unless (-e $index_file);
	};
	throw("indexing failed for $zip_vcf, $@\n") if $@;
				
	return $zip_vcf; 	
}


sub parse_seq_index {	
	my ($seq_index) = @_;
	
	my $lines = get_lines_from_file($seq_index);
	throw("Seq index file $seq_index does not exist") if ( !$lines || @$lines==0 );

	my (%desired_samples, %sample_pop);
	
	foreach my $line (@$lines) {
		next if ( $line =~ /FASTQ_FILE/i ); ### skip the headline
		my @bits = split(/\t/, $line);
		my $sample_name = $bits[9];	
		my $fastq = $bits[0];
		my $platform = $bits[12];
		my $analysis_grp = $bits[25];
		my $withdrawn = $bits[20];
		my $pop = $bits[10];
		$sample_pop{$sample_name} = $pop;	
		if ($withdrawn != 1 && 
			$platform =~ /ILLUMINA/i && 
			$analysis_grp =~ /low coverage/) {
				next if ( defined $input{pop} && $pop ne $input{pop} );
				$desired_samples{$sample_name} = 1;
		} 
	}
	return (\%desired_samples, \%sample_pop);
}



=pod

=head1 NAME

ReseqTrack/scripts/process/run_cortex.pl

=head1 SYNOPSIS

This script is the wrapper for Zam's cortex code. It is the core part of the run cortex pipeline. Another script 
cortex_process_manager.pl is used to check run status for individual jobs and creating new collections as input to kick
off down stream events.

Below lists the events in the run cortex pipeline; execept event1, all other 10 events calls run_cortex.pl. 

Pre-event: insert samples and population in input_string table

perl $ZHENG_RB_VC/scripts/process/insert_samples_into_input_string_tb.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index \
-pop LWK \
-run

Event 1: Copy fasta files and create se_ pe1, pe2 lists

perl $ZHENG_RB_VC/scripts/process/pre_run_cortex.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index \
-sample NA19449 \    #### from input_string table, with LWK_20120522_TEST as type
-store \
-save_collection \
-update \
-output_dir /nfs/1000g-work/G1K/scratch/zheng/run_cortex/input_fasta \
-run \ #### if you want to rsync the fasta


Event 2:  Create per sample graph:

perl $ZHENG_RB_VC/scripts/process/run_cortex.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection NA19449 \
-collection_type FASTQ_LIST \
-executable cortex_var_31_c1 \
-kmer_size 31 \
-mem_height 25 \
-mem_width 150 \
-para format=FASTQ \
-para max_read_len=200 \
-para quality_score_threshold=10 \
-dump_binary 1 \
-output_file_type SAMPLE_CTX \
-store

Event 3: Pool uncleaned sample graph:

perl $ZHENG_RB_VC/scripts/process/run_cortex.pl  $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection LWK_2_samples_ctx_pool \
-collection_type SAMPLE_CTX_TO_POOL \
-executable cortex_var_31_c1 \
-kmer_size 31 \
-mem_height 28 \  ### HIGH
-mem_width 100 \
-dump_covg_distribution \
-store \
-update \
-dump_binary 1 \
-output_file_type POOLED_CTX

Leave out -store when just testing construction of the command line 
-save_collection \ 
#-colour_list ./colour_list \


Event 4: cleanup the pool
## user input cleaning threshold is provided in collection name!!

perl $ZHENG_RB_VC/scripts/process/run_cortex.pl  $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection LWK_2_samples.k31.uncleaned_pool.user_input_T5 \
-collection_type POOLED_CTX \
-executable cortex_var_31_c1 \
-kmer_size 31 \
-mem_height 28 \
-mem_width 100 \
-dump_binary 1 \
-dump_covg_distribution \
-store \
-output_file_type POOL_CTX_CLEAN \
-save_collection \
-update

#--colour_list ./uncleaned_pool.colour_list \
#-para remove_low_coverage_supernodes=1 \

Event 5 :  clean up per sample graph based on cleaned pooled graph:

perl $ZHENG_RB_VC/scripts/process/run_cortex.pl  $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection NA19449.uncleaned.q10.k31.ctx \
-collection_type SAMPLE_CTX_TO_CLEAN \
-executable cortex_var_31_c2 \
-kmer_size 31 \
-mem_height 25 \
-mem_width 150 \
-para load_colours_only_where_overlap_clean_colour=0 \
-para successively_dump_cleaned_colours=CLEAN \
-output_file_type CLEANED_PER_SAMPLE_CTX \
-store

Collection SAMPLE_CTX_TO_CLEAN would contain SAMPLE_CTX and a POOL_CTX_CLEAN
#-colour_list ./cleanup_per_sample.colour_list \

Event 6: Site discovery:

perl $ZHENG_RB_VC/scripts/process/run_cortex.pl  $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection LWK_2_samples.k31.t1.cleaned_pool.ctx \
-collection_type POOL_CTX_CLEAN \
-executable cortex_var_31_c2 \
-kmer_size 31 \
-mem_height 25 \
-mem_width 150 \
-ref_binary /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37.proper_chroms.k31.ctx \
-detect_bubbles1 1/1 \
-para exclude_ref_bubbles=1 \
-para ref_colour=0 \
-output_file_type BUBBLES \
-store \
-save_collection

#-colour_list site_discovery.colour_list \


Event 7: Slice out bubble branch fasta
perl $ZHENG_RB_VC/scripts/process/make_branch_fasta.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection LWK_2_samples.k31.t1.bubbles \
-collection_type BUBBLES \
-store \
-output_file_type BUBBLE_FASTA \
-save_collection \
-kmer_size 31	


Event 8: Create bubble graph
perl $ZHENG_RB_VC/scripts/process/run_cortex.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection LWK_2_samples.k31.t1.bubbles.branches.fasta \
-collection_type BUBBLE_FASTA \
-executable cortex_var_31_c1 \
-kmer_size 31 \
-mem_height 23 \
-mem_width 100 \
-para format=FASTA \
-para max_read_len=15000 \
-dump_binary 1 \
-output_file_type BUBBLE_CTX \
-store \
-save_collection


Event 9: Overlay cleaned sample graph with bubble graph
One of the collections should be ref_binary and bubble_ctx, so this event takes care of reference bubble as well!

9a:

perl $ZHENG_RB_VC/scripts/process/run_cortex.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection NA19449.q10.k31.list_CLEAN.ctx \
-collection_type SAMPLE_CTX_CLEAN_TO_OL_BUBBLE \
-executable cortex_var_31_c2 \
-kmer_size 31 \
-mem_height 23 \
-mem_width 100 \
-para load_colours_only_where_overlap_clean_colour=0 \
-para successively_dump_cleaned_colours=OL_BUBBLES \
-output_file_type SAMPLE_CTX_CLEAN_BUBBLED \
-store

9b - ref:

perl $ZHENG_RB_VC/scripts/process/run_cortex.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection ref.ctx \
-collection_type SAMPLE_CTX_CLEAN_TO_OL_BUBBLE \
-executable cortex_var_31_c2 \
-kmer_size 31 \
-mem_height 23 \
-mem_width 100 \
-para load_colours_only_where_overlap_clean_colour=0 \
-para successively_dump_cleaned_colours=OL_BUBBLES \
-output_file_type SAMPLE_CTX_CLEAN_BUBBLED \
-store

Event 10.  Make multi-colour graph that include ref binary as well
perl $ZHENG_RB_VC/scripts/process/run_cortex.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection  LWK_2_samples.k31.t1.bubbles.ctx_pool \
-collection_type SAMPLE_CTX_CLEAN_BUBBLED_POOL \
-executable cortex_var_31_c3 \
-kmer_size 31 \
-mem_height 23 \
-mem_width 100 \
-dump_binary 1 \
-store \
-output_file_type MULTI_COLOUR_CTX \
-save_collection



Event 11: Genotyping

perl $ZHENG_RB_VC/scripts/process/run_cortex.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection  LWK_2_samples.k31.t1.bubbles.ctx_pool.multiColour_ctx \
-collection_type MULTI_COLOUR_CTX \
-executable cortex_var_31_c3 \
-kmer_size 31 \
-mem_height 23 \
-mem_width 100 \
-para max_read_len=15000 \
-para genome_size=3000000000 \
-para experiment_type=EachColourADiploidSampleExceptTheRefColour \
-para estimated_error_rate=0.01 \
-para ref_colour=0 \
-store \
-output_file_type GENOTYPE \
-save_collection


Events 13-16. Post processing, handled by other scripts

