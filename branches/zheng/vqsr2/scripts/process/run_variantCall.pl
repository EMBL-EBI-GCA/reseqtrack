#!/sw/arch/bin/perl -w 

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::HostUtils qw(get_host_object);
use ReseqTrack::Tools::GeneralUtils;
use Getopt::Long;
use ReseqTrack::Tools::BamUtils;
use File::Basename;
use ReseqTrack::Tools::Loader::File;
use Time::Local;

my $start_time = time();

$| = 1;

#### FIXME< chrom X and Y should be treated differently

my (
	%input, # config file	
	$module_name,
    $chrom,
    $region,
   	@bams,
   );
 
my $module_path = 'ReseqTrack::Tools::RunVariantCall';

GetOptions(
	\%input,    
	'dbhost=s',   
	'dbname=s',      
	'dbuser=s',
	'dbpass=s', 
	'dbport=s',  

  	'cfg_file=s',
  	       
	'collection_name=s',
	'collection_type=s',
  	'file_name_pattern=s',
  	'bam_list=s',
  	'super_pop:s', #this will query the db for chromosomal BAM files with type TD_BAM, belong to a super population (or all population), will call chrom chunk by chrom chunk
  	'bam_type:s',
  
	'program=s',
	'tabix_dir=s',
  	'reference:s',
  	'algorithm=s',
  	'parameters=s',

  	'chrom_chunk:s', #in the format of chrom:start-end, this will come from input_string table when the script is run as an event
  	'only_chrom:s', # only work on a specified chromosome, even though the input_string table contains all chromosomes; other chrs jobs would be failed
  	  	  	  	
  	'output_dir=s',
  	'output_file_type:s',
  	'output_name_prefix=s',

  	'store!',
  	'save_collection!',
 	'update!',
 	'host:s',
  	'keep_intermediate_file!',
  	'help!',
);


if ( defined $input{cfg_file} ) {
  get_params( $input{cfg_file}, \%input );
}

### Set Default
$input{update} = 0 if (!$input{update});

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd` ;
	chomp $input{output_dir};
}

if (!$input{tabix_dir}) {
	$input{tabix_dir}= "/nfs/1000g-work/G1K/work/bin/tabix/";
}	

$input{host} = '1000genomes.ebi.ac.uk' if (!$input{host});
$input{output_file_type} = "DCC_VCF" if ( !$input{output_file_type} );
$input{output_name_prefix} = "RESULTS" if ( !$input{output_name_prefix} );
$input{keep_intermediate_file} = 0 if (!$input{keep_intermediate_file});

if ($input{help}) {
	help_info();
}

print "chr chunk is $input{chrom_chunk}\n";
print "output dir is $input{output_dir}\n";
print "algorithm is $input{algorithm}\n"; 
print "super_pop is $input{super_pop}\ttype is $input{bam_type}\n";
print "parameters are $input{parameters}\n";

if ( !$input{algorithm} || $input{algorithm} !~ /samtools|umake|gatk/i ) {
	throw("Please provide snp calling algorithm you wish to use, can be samtools, gatk, or umake\n");
}

if ( !$input{reference}) {
	warn("No reference genome is provided, will use the default reference genome");
}	

if(!($input{collection_name} && $input{collection_type}) && !$input{bam_list} && !$input{file_name_pattern} && !($input{super_pop} && $input{bam_type}) ){
  throw("Please provide some input arguments by one of the following ways: " .
  		"-collection_name and -collection_type\n" .
  		"-bam_list\n" .
  		"-file_name_pattern\n" .
        "-super_pop and -bam_type\n" );
}

### Fail jobs with chrom_chunk from input_string table that don't match $only_chrom 
if ($input{chrom_chunk}) {
	if ($input{chrom_chunk} =~ /:/ ) {
		($chrom, $region) = split(/:/, $input{chrom_chunk});
	}
	else {
		$chrom = $input{chrom_chunk};
	}	
	throw("chrom_chunk $input{chrom_chunk} does not belong to chromosome specified by -only_chrom $input{only_chrom}") if ($input{only_chrom} && ($input{chrom} ne $input{only_chrom})); ## "ne" works for both strings and numbers
}	

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
    );

my $fa = $db->get_FileAdaptor;

if($input{collection_type} && $input{collection_name}){
  my $collection =
    $db->get_CollectionAdaptor->fetch_by_name_and_type($input{collection_name},
                                                        $input{collection_type});
  throw("No collection found for $input{collection_name} and $input{collection_type}") if (!$collection); 
  foreach my $file (@{$collection->others}){
    next if ( check_file_sanity($file->name, $chrom) eq 'skip' );
    push(@bams, $file->name);
  }
}

if($input{bam_list}){
  my $lines = get_lines_from_file($input{bam_list});
  foreach my $line(@$lines){
    next if ( check_file_sanity($line, $chrom) eq 'skip');
    push(@bams, $line);
  }
}

if($input{file_name_pattern}){
  my $files = $fa->fetch_all_like_name($input{file_name_pattern});
  throw("No files found containing pattern $input{file_name_pattern}") if (!$files || @$files ==0);
  foreach my $file(@$files){
    next if ( check_file_sanity($file->name, $chrom) eq 'skip' );
    push(@bams, $file->name);
  }
}

my $super_pop_name = "unknownSupPop"; ## set default
if ($input{super_pop} && $input{bam_type}) {
	my $files = $fa->fetch_by_type($input{bam_type});
	throw("No files found with type $input{bam_type}") if (!$files || @$files==0);
	if ($input{super_pop} =~ /all/i) {
		$super_pop_name = "all";
		foreach my $file(@$files) {
			next if ( check_file_sanity($file->name, $chrom) eq 'skip' );
 			push(@bams, $file->name);
  		}
	}
	else {  
		throw("super population input should be 'all' or super_pop_name:string, string is sub pop seperated by ','") if ($input{super_pop} !~ /,|:/); 
		my $string;
		($super_pop_name, $string) = split(/:/, $input{super_pop});
		my @pops = split(/,/, $string);  #super_pop is a string, led by super_pop name, followed by sub-pops separated by ','  EUR:IBS,CLM,CEU
		foreach my $file(@$files) {
			next if ( check_file_sanity($file->name, $chrom) eq 'skip' );
			foreach my $pop ( @pops ) {
				if ( $file->name =~ /$pop/i) {
				 	push(@bams, $file->name);
				}	
			}	
  		}
	}	
}

$db->dbc->disconnect_when_inactive(1);

print "Input BAMs are:\n";
print join ("\n", @bams), "\n";

if ($input{algorithm} =~ /samtools/i) {
	$module_name = "CallBySamtools";
}	
elsif ( $input{algorithm} =~ /gatk/i || $input{algorithm} =~ /unified_genotyper/i ) {
	$module_name = "CallByGATK";
}
elsif ( $input{algorithm} =~ /umake/i ) {
	$module_name = "CallByUmake";
}	
else {
	throw("Please use one of the 3 algorithms to make variant calls - samtools, gatk or umake\n");
}

my $varCalling_module = $module_path . '::'. $module_name;

my $constructor_hash;
$constructor_hash->{-parameters} = parameters_hash($input{parameters});
$constructor_hash->{-input_files} = \@bams;
$constructor_hash->{-program} = $input{program};
$constructor_hash->{-working_dir} = $input{output_dir};
$constructor_hash->{-reference} = $input{reference};
$constructor_hash->{-chrom} = $chrom if ($chrom);
$constructor_hash->{-region} = $region if ($region);
$constructor_hash->{-save_files_for_deletion} = $input{keep_intermediate_file};
$constructor_hash->{-output_name_prefix} = $input{output_name_prefix};
$constructor_hash->{-super_pop_name} = $super_pop_name; 

my $object = construct_new_object($varCalling_module, $constructor_hash);

$object->run;

my $flt_vcf = $object->output_files;

$db->dbc->disconnect_when_inactive(0);

foreach my $outfile ( @$flt_vcf ) {  
	if ( -e $outfile ) {	
		my $loader;
		if ($outfile =~ /vcf$|vcf.gz$/ ) {
			my $zipped_vcf = bgzip_and_index($outfile, $input{tabix_dir});
			$outfile = $zipped_vcf;
			my $tabix_file = $outfile . ".tbi";
		
			$loader = ReseqTrack::Tools::Loader::File->new
		 		 (
				   -file => [$outfile, $tabix_file],
				   -do_md5 => 1,
				   -hostname => $input{host},
				   -db => $db,
				   -assign_types => 0,
				   -check_types => 0,
				   -type => $input{output_file_type},
				   -update_existing => $input{update},
				  );
		}

		next if ($outfile =~ /vcf.gz.tbi$/);
		
		my $file_objects;
		if($input{store}){
		  $loader->process_input();
		  $loader->create_objects();
		  $loader->sanity_check_objects();
		  $file_objects = $loader->load_objects();
		}
		
		### save the VCF file in collection_by_chrom, 
		### after many jobs of the event are run, each collection would contain vcf files of different chrom chunks
		if ($input{save_collection} && $outfile !~ /tbi$/i ) {
			my $collection_by_chrom_name = $input{output_name_prefix} . "_" . $input{algorithm} . "_$super_pop_name" . "_chr$chrom";
			if ($outfile =~ /vcf/) {
				my $collection_by_chrom =  ReseqTrack::Collection->new(
				  -name => $collection_by_chrom_name,
				  -others => $file_objects,
				  -type => $input{output_file_type},  
				  -table_name => 'file',
				);
				$db->get_CollectionAdaptor->store($collection_by_chrom);  ## new files will be added in without over write the old ones.
			}
		}	
		
		print "*********************************************************************************************************************************\n";
		print "**** Output filtered VCF file path is $outfile\n";
		if ($input{store}) {
			print "**** Output filtered VCF file $outfile has been stored in the database\n";
		}
		else {
			print "**** Output filtered VCF file path is $outfile; it is NOT stored in the database as -store is not specified\n";
		}	
		print "*********************************************************************************************************************************\n";  
	}	
	else {
		throw("Output file " . $outfile . " does not exist");
	}	
}

my $end_time = time();
print "Job took " . ($end_time - $start_time)/60 . " minutes\n";

##############
#### SUBS ####
##############
sub parameters_hash{
  my ($string) = @_;

  my %parameters_hash;

  if ($string) {
      
    my @pairs;
    if ($string !~ /,/ && $string =~ /:/ ) {
		push @pairs, $string;
    }       
    elsif($string =~  /,/ && $string =~ /:/){
		@pairs = split (/,/, $string);
    }          
    else {
          throw("Please provide running parameters as name and value separated by ':'  Multiple name and value pairs can be used and should be separated by ',' ");
	}
	
	foreach my $pair(@pairs){
		my ($key, $value) = split (/:/, $pair);
			$key   =~ s/^\s+|\s+$//g;
          	$value =~ s/^\s+|\s+$//g; 
          	$parameters_hash{$key} = $value;
          	print "key is $key, value is $value\n";	 
    }#end of foreach pairs
  }# end of if ($string)
  
  return \%parameters_hash;
}  		 		
  		
sub construct_new_object {
  my ($module, $args) = @_;
  my $file = "$module.pm";
  $file =~ s{::}{/}g;
  eval {
    require "$file"; # raises exception upon error, checks @INC
  };
  if($@){
    throw("ReseqTrack::Tools::EventPipeline::setup_batch_submission_system ".
          "Can't find $file [$@]");
  }
  my %constructor_args = %$args;
  my $object = $module->new(%constructor_args);
  return $object;
}

sub check_file_sanity{
  my ($file, $chrom) = @_;
  my $flag = 'ok';
 
  throw($file." does not exist on the current filesystem") unless(-e $file);
  $flag = 'skip' if ($file !~ /\.bam$/);
 
  if ($chrom) {
      
      my @bits = split(/_|\./, $file);
      my $file_chr = "NA";
      foreach my $bit ( @bits) {
          if ($bit =~ /chr/i) {
              $file_chr = $bit;
              $file_chr =~ s/chrom//i;
              $file_chr =~ s/chr//i;
          }
      }        

      $flag = 'skip' if ( $file_chr ne $chrom );
  }     
  return $flag;
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
	
			    
sub help_info {

 exec('perldoc', $0);

}


#############################################################################################

=pod

=head1 NAME

	~/ReseqTrack/scripts/process/run_variantCall.pl

=head1 Aarguments:

	cfg_file,			This specifies a configuration file that contains different arguments to pass to the program as listed below:

 Database arguments (required):
	
	dbhost, 			the name of the mysqlhost
	dbname, 			the name of the mysql database
	dbuser, 			the name of the mysql user
	dbpass, 			the database password if appropriate
	dbport, 			the port the mysql instance is running on, this defaults to 4197 the standard
						port for mysqlg1kdcc.ebi.ac.uk
    
 Run program arguments (required):
    
	algorithm,			can be one of the three: samtools, gatk, umake
	reference,			Reference genome; default is ncbi build 37  (umake requires a umfa file in addition to the fa file)		
	
 BAM input arguments (4 ways to input BAMs for calling, one and only one is required):
	
	a.		 
	collection_name,		BAMs belong to a BAM collection with this name and this type
	collection_type,		Type of BAM collection
	
	b.
	file_name_pattern,		BAMs with this text string in their names
	
	c.
	bam_list,			BAMs in this file that lists full BAM pathes
	
	d.
	super_pop,			This option is often used for running the script as farm jobs.  
						DB will be queried to find transformed chromosomal BAM files (type TD_BAM) that belong to a super population (or all populations)
						If call for all population, use "all";
						If call for a specific super population, use this format:
						super pop name:pop1,pop2,pop3 (example EUR:IBS,CLM,CEU)
	bam_type,			type needs to be one of the transformed BAMs such as TD_BAM or TD_PHASE1_BAM
	
 Chromosome or chromosomal region arguments (optional):
	
	chrom_chunk, 			in the format of chrom:start-end, this will come from input_string table when the script is run as an event
						it can also specify an entire chromosome by simply leaving out ":start-end".
						optional for samtools and gatk but mandatory for umake (both chrom and region are mandatory)
	only_chrom,  			only work on a specified chromosome, even though the input_string table contains all chromosomes; other chrs jobs would be failed
	
 Output arguments (optional):
	
	output_dir,			This is where all temporary files and output files will be written to
	output_name_prefix,		A short word you may like to put in the output vcf file and collection name as prefix to identify this round of calling 
						(example "phase1" or "mouse" or "Tony")
	store,				To store the output VCF file in the database
	save_collection,		To store the output file into appropriate collections 
	update,				To replace an output VCF file in the database
	keep_intermediate_file		Do not delete intermediate files
	
 Variant calling parameters (required):
	
	parameters,			Key and value pairs to pass options to different variant calling algorithms. Options are different for different callers.  Key and value are separated by ':'; if there are multiple key and value pairs, they 
should be separated by ','  The pipeline does not check for the accuracy of the input parameter keys so take extra caution to use the precise string as defined by different callers.  
Below lists the options we have implemented for each caller to take now. Please let us know if you want additional input parameters to be implemented.

		Samtools:
		
		Commonly used parameter keys can be one or more of the followings. Please see Samtools page for details -
		http://samtools.sourceforge.net/samtools.shtml#3
		
		mpileup:	string, command line options to use with "samtools mpileup"
      				Here is a list of options: samtools mpileup [-EBug] [-C capQcoef]  [-l list] [-M capMapQ] [-Q minBaseQ] [-q minMapQ]
		bcfview: 	string, command line options to use with "bcftools view"
      				Here is a list of options: bcftools view [-AbFGNQSucgv] [-D seqDict] [-l listLoci] [-i gapSNPratio] [-t mutRate] [-p varThres] [-P prior] [-1 nGroup1] [-d minFrac] [-U nPerm] [-X permThres] [-T trioType]
		vcfutils:	string, command line options to use with "vcfutils.pl"
		bcftools_path: 
		vcfutils_path:
		
		For mpileup:, values can be taken from:
		samtools mpileup [-EBug] [-C capQcoef] [-l list] [-M capMapQ] [-Q minBaseQ] [-q minMapQ]
		The [-r region] and [-f reference] have been coded as separate input entries so don't include them in the parameteres
		
		For bcfview:, values can be taken from:
		bcftools view [-AbFGNQSucgv] [-D seqDict] [-l listLoci] [-i gapSNPratio] [-t mutRate] [-p varThres] [-P prior] [-1 nGroup1] [-d minFrac] [-U nPerm] [-X permThres] [-T trioType]  
		[-s listSample] and [region] have been coded separately so don't include them in the parameters
		
		For vcfUtils:, 
		The -D option of varFilter controls the maximum read depth, which should be adjusted to about twice the average read depth. 
		
		gatk:
		
		Commonly used parameter keys can be one or more of the followings. Please see GATK website for detailed descriptions and additional parameters - 
		http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_genotyper_UnifiedGenotyper.html
		dbSNP: 
  		dcov: 				integer, command line options to use with "samtools mpileup"
  		computeSLOD:		If provided, we will calculate the SLOD. This argument is not enabled by default because it increases the runtime by an appreciable amount.
  		stand_emit_conf:	float, The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be emitted (and filtered if less than the calling threshold).
  		stand_call_conf:	float, The minimum phred-scaled confidence threshold at which variants not at 'trigger' track sites should be called. 
  		debug_file:			File to print all of the annotated and detailed debugging output.
  		glm:				string, genotype_likelihood_model, choose one of the followings "SNP", "BOTH", "INDEL"
  		genotyping_mode:	String, Default is DISCOVERY. Should we output confident genotypes (i.e. including ref calls) or just the variants? Choose one of the two "DISCOVERY" or "GENOTYPE_GIVEN_ALLELES"
  		pcr_error_rate:		Double with default value 1.0E-4
  		p_nonref_model:		Non-reference probability calculation model to employ -- EXACT is the default option, while GRID_SEARCH is also available. 
  		minIndelCnt:		Minimum number of consensus indels required to trigger genotyping run. int with default value 5
  		output_mode:		Should we output confident genotypes (i.e. including ref calls) or just the variants? Default is EMIT_VARIANTS_ONLY
  		min_mapping_quality_score:	Minimum read mapping quality required to consider a read for calling.
  		min_base_quality_score:		Minimum base quality required to consider a base for calling.
  		metrics_file:		File to print any relevant callability metrics output.
  		max_deletion_fraction:		Maximum fraction of reads with deletions spanning this locus for it to be callable [to disable, set to < 0 or > 1; default:0.05].
  		indel_heterozygosity:		Heterozygosity for indel calling. This argument informs the prior probability of having an indel at a site.
  		heterozygosity:		Heterozygosity value used to compute prior likelihoods for any locus.
  		group:				One or more classes/groups of annotations to apply to variant calls. Which groups of annotations to add to the output VCF file.
  		annotation:
  		alleles:			The set of alleles at which to genotype when in GENOTYPE_MODE = GENOTYPE_GIVEN_ALLELES.
		
		umake:

		Commonly used parameters keys can be one or more of the following. Additional parameters can be found in UMAKE website - 
		http://genome.sph.umich.edu/wiki/UMAKE#Configuration_File
				
		offset_off_target:	option for EXOME sequencing, extend target by given # of bases  
		target_bed:  
  		dbSNP: 			string, path to prefix of dbsnp vcf rod file
  		indel_prefix:	string, path to prefix of indel file
  		hm3_prefix: 	string, path to prefix of hapmap3 file
  		FILTER_MAX_SAMPLE_DP:	argument for filtering
  		FILTER_MIN_SAMPLE_DP:	argument for filtering
  		FILTER_MQ:		filter output calls by samtools view MQ, default is 3
  		FILTER_FLAG:	filter output calls by samtools view flag, deafult is 0x0704
  		unit_chunk:		size of chromosomal chunk to call variants on; default is 5Mb
  		LD_Nsnps:		Chunk size of genotype refinement, default is 10,000 SNPs
  		LD_overlap:		overlapping sizes of chunks, default is 1,000 SNPs
		

=head1 SYNOPSIS

	This is a generic script to call variants from BAM files. Algorthms that can be used for making varaint call include: 
	samtools (http://samtools.sourceforge.net/mpileup.shtml)   
	gatk (http://www.broadinstitute.org/gsa/wiki/index.php/Unified_genotyper)
	umake (http://genome.sph.umich.edu/wiki/UMAKE)
	
	Input BAMs are provided by 4 different ways as described above. Output is variants discovered in the format of VCF files. 

=head1	OUTPUT

	Output files from samtools and gatk can be found in the user-specified output_dir; if an output_dir is not specified, the current 
	directory is the output_dir.
	
	Output VCF files are organized in sub-directories of different chromosomes. 
	
	Output file name contains the following elements:
	a user-specified prefix (such as PHASE1, TEST, TONY), default is RESULTS 
	a super population name such as EUR or ALL taken from -super_pop option, if this option not available, default is unknownSupPop;
	number of BAMs used for calling
	chromosomal region if specified
	algorithm name
	
	Below is a few examples of VCF files:
	TEST_unknownSupPop_of_1bams.chr10_1000000-2000000.samtools.flt.vcf.gz
	RESULTS_unknownSupPop_of_1bams.chr10_10000000-10500000.gatk.vcf.gz
	
	Output for UMAKE is defined by UMAKE itself; it contains lots of intermediate files organized in diffrent sub-directories specified 
	in section "RELATIVE DIRECTORY UNDER OUT_DIR" of the configuration file and separated into different chromosomes. For --snpcall, the 
	most useful output is chrN.filtered.vcf.gz file under vcfs/chrN; to distinguish vcf files generated by different chromosome chunks,
	we rename the file to include chromosome chunk info, such as in chr20.1000000-2000000.filtered.vcf.gz. When keep_intermediate_file is et to '0', 
	all but the filtered.vcf.gz will be kept. 
	
	The script run_variantCall.pl will print the path to the most relevant output file when a run is done.  Output VCF files will be stored 
	in the database as type "DCC_VCF" if option '-store' is specified. They will also be stored in specific collections of type "DCC_VCF" 
	if -save_collection is specified.

=head1 EXAMPLE COMMAND LINE

perl $ZHENG_RT/scripts/process/run_variantCall.pl \
-cfg_file /nfs/1000g-work/G1K/work/zheng/snp_calling/samtools.cfg -chrom_chunk 20:5000000-5500000 -save_collection


=head1 EXAMPLE SAMTOOLS CONFIG FILE


 ########## Calling algorithm and parameters ##########
 algorithm=samtools
 reference=/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta
 parameters=mpileup:-EDS -C50 -m2 -F0.0005 -d 10000 -P ILLUMINA -ug,bcfview:-p 0.99 -bvcg,vcfutils:-D 10000 -d 2 -a 2 -Q 10 -w 10 -W 3 -1 1e-5 -4 1e-7

 ########## Input BAM files ##########
 #super_pop=all
 #super_pop=TEST:ASW,
 #bam_type=TD_PHASE1_BAM
 #only_chrom=20

 OR

 # Test case 1 (a quick test case)
 collection_name=NA06985.ILLUMINA.bwa.low_coverage
 collection_type=TEST_BAM
 chrom_chunk=10:1000000-1050000

 OR

 # Test case 2
 file_name_pattern=NA06985.chrom11.ILLUMINA.bwa.CEU.low_coverage.2011b

 OR

 # Test case 3
 bam_list=/nfs/1000g-work/G1K/work/zheng/snp_calling/bam_list

 ########## Output dir ##########
 output_dir=/nfs/1000g-work/G1K/work/zheng/snp_calling/samtools/
 #output_name_prefix='TEST'
 keep_intermediate_file=1
 store=1
 update=1
 save_collection=1

 ########## database ##########
 dbhost=mysql-g1kdcc-public
 dbuser=g1krw
 dbpass=thousandgenomes
 dbport=4197
 #dbname=g1k_archive_staging_track
 dbname=zheng_var_call

=head1 EXAMPLE GATK CONFIG FILE

 ########## Calling algorithm and parameters ##########
 algorithm=gatk
 reference=/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/human_g1k_v37.fasta
 parameters=dcov:40,stand_emit_conf:2.0,stand_call_conf:4.0,glm:BOTH,dbSNP:/nfs/1000g-work/G1K/scratch/zheng/reference_genomes/dbsnp132_20101103.vcf.gz

 ########## Input BAM files ##########
 #test case1:
 bam_list=/nfs/1000g-work/G1K/work/zheng/snp_calling/bam_list
 chrom_chunk=10:10000000-10500000

 OR

 #test case2 (a quick test):
 #collection_name=NA06985.ILLUMINA.bwa.low_coverage
 #collection_type=TEST_BAM
 #chrom_chunk=10:1000000-2000000

 OR

 #run as event case:
 #super_pop=all
 #super_pop=TEST:ASW,
 #bam_type=TD_PHASE1_BAM
 #only_chrom=20

 ########## Output dir ##########
 output_dir=/nfs/1000g-work/G1K/work/zheng/snp_calling/gatk/results
 keep_intermediate_file=1
 store=1
 update=1
 save_collection=1

 ########## database ##########
 dbhost=mysql-g1kdcc-public
 dbuser=g1krw
 dbpass=thousandgenomes
 dbport=4197
 #dbname=g1k_archive_staging_track
 dbname=zheng_var_call


=head1 EXAMPLE UMAKE CONFIG FILE


 ########## Calling algorithm and parameters ##########
 algorithm=umake
 reference=/nfs/1000g-work/G1K/work/bin/umake-resources/ref/human.g1k.v37.fa
 parameters=FILTER_MQ:3,dbSNP:/nfs/1000g-work/G1K/work/bin/umake-resources/dbSNP/dbsnp_129_b37.rod,hm3_prefix:/nfs/1000g-work/G1K/work/bin/umake-resources/HapMap3/hapmap3_r3_b37_fwd.consensus.qc.poly,indel_prefix:/nfs/1000g-work/G1K/work/bin/umake-resources/indels/1kg.pilot_release.merged.indels.sites.hg19

 ########## Input BAM files ##########
 # Test case 1
 bam_list=/nfs/1000g-work/G1K/work/zheng/snp_calling/bam_list
 chrom_chunk=10:10000000-13000000
 
 OR 

 # run as event
 #super_pop=all
 #super_pop=TEST:ASW,
 #bam_type=TD_PHASE1_BAM
 #only_chrom=20

 ########## Output dir ##########
 output_dir=/nfs/1000g-work/G1K/work/zheng/snp_calling/umake/results
 keep_intermediate_file=1
 store=1
 save_collection=1
 update=1

 ########## database ##########
 dbhost=mysql-g1kdcc-public
 dbuser=g1krw
 dbpass=thousandgenomes
 dbport=4197
 #dbname=g1k_archive_staging_track
 dbname=zheng_var_call


 
 
