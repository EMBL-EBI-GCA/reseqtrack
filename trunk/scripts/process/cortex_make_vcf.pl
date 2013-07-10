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
	'log_dir=s',
	'sample_name_list_dir=s',
	'ref_genome=s',
	
  	'output_dir=s',
  	'output_file_type:s',
  	'tabix_dir:s',

  	'store!',
  	'save_collection!',
 	'update!',
 	'host:s',
  	'run!',		
);

if ( !$input{collection} ) {
	throw("Please provide a collection name\n");
}
if ( !$input{log_dir} ) {
	throw("Please provide a log dir where genotype log was written to\n");
}

if (!$input{tabix_dir}) {
	$input{tabix_dir}= "/nfs/1000g-work/G1K/work/bin/tabix/";
}	

$input{host} = '1000genomes.ebi.ac.uk' if (!$input{host});
$input{update} = 0 if (!$input{update});
$input{collection_type} = "SAM" if (!$input{collection_type});
$input{output_file_type} = "VCF" if (!$input{output_file_type});
$input{sample_name_list_dir} = "/nfs/1000g-work/G1K/work/zheng/cortex/post_process/" if ( !$input{sample_name_list_dir} );
$input{ref_genome} = "/nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37.fasta.proper_chroms_only" if (!$input{ref_genome});

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

my $col = $db->get_CollectionAdaptor->fetch_by_name_and_type($input{collection}, $input{collection_type});
my $sam_files = $col->others;

$db->dbc->disconnect_when_inactive(1); 

foreach my $sam ( @$sam_files ) {
	throw("sam file " . $sam->name . " does not exist") unless (-e $sam->name);
	
	my $outd = $input{output_dir} . "/$pop";
	mkpath($outd) unless (-e $outd);
	
	my @tmp = split(/\_/, basename($sam->name));
	my $pop = $tmp[0];
	my $num_colours = $tmp[1] + 1;

	my $genotype_file = query_file_by_type_and_pop("GENOTYPE", $pop);
	throw("Genotype file does not exist") unless (-e $genotype_file);
	
	my $genotype_log = get_geno_log($input{log_dir});
	
	my $sample_name_list = $input{sample_name_list_dir} . $pop . ".sample_name_list"; 
	throw("sample name list $sample_name_list does not exist") unless (-e $sample_name_list);
	
	my $classified_var_file = query_file_by_type_and_pop("CLASSIFIED", $pop);
	throw("Classified variant file does not exist") unless (-e $classified_var_file);
	
	my $out_vcf_prefix = $pop . "_" . $tmp[1] . "_samples";
	my $output_decomp_vcf = $outd . "/" . $out_vcf_prefix . ".decomp.vcf";
	my $output_raw_vcf = $outd . "/" . $out_vcf_prefix . ".raw.vcf";
	
	throw("Output vcf file $output_decomp_vcf and/or $output_raw_vcf exist") if (-e $output_decomp_vcf || -e $output_raw_vcf);
	
	my $command = "perl  /nfs/1000g-work/G1K/work/bin/cortex/scripts/analyse_variants/process_calls.pl ";
	$command .= "--callfile $genotype_file ";
	$command .= "--callfile_log $genotype_log ";
	$command .= "--outvcf $out_vcf_prefix ";
	$command .= "--outdir  $outd ";
	$command .= "--samplename_list $sample_name_list ";
	$command .= "--num_cols $num_colours ";
	$command .= "--stampy_hash /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37 ";
	$command .= "--vcftools_dir /nfs/1000g-work/G1K/work/bin/vcftools/ ";
	$command .= "--caller BC ";
	$command .= "--kmer 31 ";
	$command .= "--stampy_bin /nfs/1000g-work/G1K/work/bin/stampy/stampy.py "; ## This may be left out
	$command .= "--refcol 0 ";
	$command .= "--pop_classifier $classified_var_file ";
	$command .= "--ploidy 2 ";
	$command .= "--ref_fasta $input{ref_genome} ";
	
	print "Run ...\n$command\n";
	`$command` if $input{run};	
		
	if (-e $output_decomp_vcf && -e $output_raw_vcf ) {
		my $zipped_decomp_vcf = bgzip_and_index($output_decomp_vcf, $input{tabix_dir});
		my $zipped_raw_vcf = bgzip_and_index($output_raw_vcf, $input{tabix_dir});		
		save_file($zipped_decomp_vcf, "VCF");
		save_file($zipped_raw_vcf, "VCF");
	}
	else {
		throw("Output VCF files $output_decomp_vcf and $output_raw_vcf do not exist\n");
	}		

}	

############ SUBS ##########
sub get_geno_log {
	my ($log_dir) = @_;
	
	my $genotype_logs = `find $log_dir -name \"\*$pop*genotype\*out\" -print`;
	my @logs = split(/\n/, $genotype_logs);
	my $cnt = 0;
	my $genotype_log;
	foreach my $log ( @logs ) {
		next if $log !~ /event11/;
		print "Geno log file is\n$log\n";
		$cnt++;
		$genotype_log = $log;
	}	
	throw("More than one genotype log files found for $input{pop} in $input{log_dir}") if $cnt > 1;
	throw("No genotype log file found for $input{pop} in $input{log_dir}") if ($cnt == 0);
	return $genotype_log;
}
	
sub query_file_by_type_and_pop {
	my ($type, $population) = @_;
	my $files = $db->get_FileAdaptor->fetch_by_type($type);
	if (!$files || @$files == 0 ) {
		throw("No file objs found for $type");
	}	
	my $cnt = 0;
	my $correct_file;
	foreach my $file ( @$files ) {
		if ( $file->name =~ /$population/ ) {
			throw($file->name . "does not exist") unless ( -e $file->name);
			$correct_file = $file->name;
			$cnt++;
		}
	}
	
	if ($cnt == 0 ) {
		throw("No file object found for $type of pop $population");	
	}	
	elsif ($cnt > 1) {
		throw("More than one files found for $type of pop $population");	
	}
	return $correct_file;	
}

sub save_collection { 
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
	
=pod

=head1 NAME

ReseqTrack/scripts/process/cortex_make_vcf.pl

=head1 SYNOPSIS

This script is the last step of run cortex pipeline. It takes the original genotype file and the classified (filtered) genotype file and 
5prime_falnking alignment files (sam), and generate the final VCF files.

=head1 EXAMPLES
  
perl $ZHENG_RB_VC/scripts/process/cortex_make_vcf.pl $WRITE_DB_ARGS -dbname zheng_run_cortex_test\
-output_dir /nfs/1000g-work/G1K/work/zheng/cortex/post_process \
-collection_type SAM \
-output_file_type VCF \
-collection PEL_50_bubbled_sample_ctxs_total2.sam \
-log_dir /nfs/nobackup/resequencing_informatics/zheng/run_cortex/log






__DATA__
=head	
bsub -R "rusage[mem=80000]" -M80000 -q 1000genomes -o /nfs/1000g-work/G1K/work/zheng/cortex/post_process/GIH_78_samples.convert_to_vcf.log2 \

perl  /nfs/1000g-work/G1K/work/bin/cortex/scripts/analyse_variants/process_calls.pl \
--callfile /nfs/1000g-work/G1K/archive_staging/ftp/technical/working/20121015_cortex_phase2/GIH/GIH_bubbles_thresh4.genotyped.20121022 \
--callfile_log   /nfs/nobackup/resequencing_informatics/zheng/run_cortex/log/event11/78/GIH_78_bubbled_sample_ctxs.multiColour_ctx_genotype_7422.out \
--outvcf GIH_78_samples.vcf \
--outdir  /nfs/1000g-work/G1K/work/zheng/cortex/post_process/ \
--samplename_list /nfs/1000g-work/G1K/work/zheng/cortex/post_process/GIH.samplename_list \
--num_cols 79 \
--stampy_hash /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37 \
--vcftools_dir /nfs/1000g-work/G1K/work/bin/vcftools/ \
--caller BC \
--kmer 31\
--stampy_bin /nfs/1000g-work/G1K/work/bin/stampy/stampy.py \
--refcol 0 \
--pop_classifier /nfs/1000g-work/G1K/work/zheng/cortex/post_process/GIH_78_bubbled_sample_ctxs.multiColour_ctx_genotype.classified \
--ploidy 2 \
--ref_fasta /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37.fasta.proper_chroms_only 


perl  /nfs/1000g-work/G1K/work/bin/cortex/scripts/analyse_variants/process_calls.pl \
--callfile /nfs/nobackup/resequencing_informatics/zheng/run_cortex/results/ACB_64_bubbled_sample_ctxs.genotype \
--callfile_log /nfs/nobackup/resequencing_informatics/zheng/run_cortex/log/event11/23/ACB_64_bubbled_sample_ctxs.multiColour_ctx_genotype_8111.out \
--outvcf ACB_64_samples \
--outdir  /nfs/1000g-work/G1K/work/zheng/cortex/post_process/ACB \
--samplename_list /nfs/1000g-work/G1K/work/zheng/cortex/post_process/ACB.sample_name_list \
--num_cols 65 \
--stampy_hash /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37 \
--vcftools_dir /nfs/1000g-work/G1K/work/bin/vcftools/ \
--caller BC --kmer 31 --stampy_bin /nfs/1000g-work/G1K/work/bin/stampy/stampy.py --refcol 0 \
--pop_classifier /nfs/1000g-work/G1K/work/zheng/cortex/post_process/ACB/ACB_8274411_65.classified \
--ploidy 2 \
--ref_fasta /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37.fasta.proper_chroms_only 


perl  /nfs/1000g-work/G1K/work/bin/cortex/scripts/analyse_variants/process_calls.pl \
--callfile /nfs/nobackup/resequencing_informatics/zheng/run_cortex/results/PEL_50_bubbled_sample_ctxs.genotype \
--callfile_log /nfs/nobackup/resequencing_informatics/zheng/run_cortex/log/event11/21/PEL_50_bubbled_sample_ctxs.multiColour_ctx_genotype_1624.out \
--outvcf PEL_50_samples \
--outdir  /nfs/1000g-work/G1K/work/zheng/cortex/post_process/PEL \
--samplename_list /nfs/1000g-work/G1K/work/zheng/cortex/post_process/PEL.sample_name_list \
--num_cols 51 \
--stampy_hash /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37 \
--vcftools_dir /nfs/1000g-work/G1K/work/bin/vcftools/ \
--caller BC \
--kmer 31 \
--stampy_bin /nfs/1000g-work/G1K/work/bin/stampy/stampy.py \
--refcol 0 \
--pop_classifier 1 \
--ploidy 2 \
--ref_fasta /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37.fasta.proper_chroms_only 
Output VCF files /nfs/1000g-work/G1K/work/zheng/cortex/post_process/PEL/PEL_50_samples.decomp.vcf.gz and /nfs/1000g-work/G1K/work/zheng/cortex/post_process/PEL/PEL_50_samples.raw.vcf.gz do not exist


=cut	
