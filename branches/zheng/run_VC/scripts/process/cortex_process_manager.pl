#!/sw/arch/bin/perl -w

use strict;
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
use ReseqTrack::Tools::Intersection;
use ReseqTrack::Tools::RunSplit;

my $start_time = time();

$| = 1;

my (
	%input,
	%samples_have_ctx,
	@uncleaned_sample_ctxs,
);		


GetOptions(
	\%input,    
	'dbhost=s',   
	'dbname=s',      
	'dbuser=s',
	'dbpass=s', 
	'dbport=s',

	'seq_index=s',
	'pop=s',	
	'cleaning_threshold:i',
	'ref_binary:s',
	'log_dir:s',
	'out_dir:s',
	'run!',

  	'save_collection!',
  	'store!',
 	'update!',
 	'host:s',
  	'help!',	
  	'post_processing!',	
);

$input{update} = 0 if (!$input{update});
$input{host} = '1000genomes.ebi.ac.uk' if (!$input{host});
if ( !$input{seq_index}  || !$input{pop}) {
	throw("Please provide a sequence index file and a population name to be processed");
}

if ( $input{post_processing} && ( !$input{log_dir} || !$input{out_dir} )) {
	throw("Please provide out_dir to where output files are written and log_dir from where genotype log file will be found"); 
}

my $out_dir_pop;
if ( $input{post_processing} ) {
	$out_dir_pop = $input{out_dir} . "/$input{pop}";
	mkpath($out_dir_pop) unless (-e $out_dir_pop);
}

if ( $input{log_dir}  && ($input{out_dir} =~ /rerun/i && $input{log_dir} !~ /rerun/i) ) {
	throw("Check if you have used the right outdir and log dir in the command line");
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
);

my ($samples_to_look_for, $sample_to_pop_mapping) = parse_seq_index($input{seq_index});
print "For pop $input{pop}, total number of samples input to cortex is " . keys(%$samples_to_look_for) . "\n";

my $total_sample_cnt = keys(%$samples_to_look_for);
#my $total_sample_cnt = 2; 
### FIXME, change this after testing

if ($input{post_processing}) {
	goto SKIP6;
}	

############################################################
### Check if event 2 is finished for all samples
### Create collection SAMPLE_CTX_TO_POOL
#############################################################
my $sample_ctxs = $db->get_FileAdaptor->fetch_by_type("SAMPLE_CTX");
foreach my $sample_ctx ( @$sample_ctxs ) {
	next if ( am_i_the_right_entry($sample_ctx) == 0);
	my @tmp = split(/\./, basename($sample_ctx->name));
	my $sam = $tmp[0];
	$samples_have_ctx{$sam} = 1 ;
	push @uncleaned_sample_ctxs, $sample_ctx;
}

my $exist_sample_ctx_to_pools = $db->get_CollectionAdaptor->fetch_by_type("SAMPLE_CTX_TO_POOL");
foreach my $exist_sample_ctx_to_pool ( @$exist_sample_ctx_to_pools ) {
	if ( am_i_the_right_entry($exist_sample_ctx_to_pool) == 1 ) {
		print "Event 2 completed, skip\n";
		goto SKIP1;
	}		
}
	
if ( keys(%$samples_to_look_for) == keys(%samples_have_ctx) ) {
	### FIXME, change != to == after testing!!!
	print "Event 2 (create per sample ctx) completed for pop $input{pop}\n";
	my $collection_name = $input{pop} . "_" . keys(%$samples_to_look_for) . "_sample_ctxs_to_pool";
	my $collection =  ReseqTrack::Collection->new(
		-name => $collection_name,
		-others => \@uncleaned_sample_ctxs,
		-type => "SAMPLE_CTX_TO_POOL",  
		-table_name => 'file',
	);
	save_collection($collection, $collection_name, "SAMPLE_CTX_TO_POOL");			
}	
else {
	print "Event 2 (create per sample ctx) has not completed for pop $input{pop} yet:\n";
	print "Number of samples to process is " .  keys(%$samples_to_look_for) . "\n";
	print "Number completed so far is " . keys(%samples_have_ctx) . "\n";
	my @all_samples = keys(%$samples_to_look_for);
	my @samples_processed = keys(%samples_have_ctx);
	my $all_samples_object = ReseqTrack::Tools::Intersection
                     ->new(
                           -LIST => \@all_samples,
                          );
    my $samples_processed_obj =  ReseqTrack::Tools::Intersection
                     ->new(
                           -LIST => \@samples_processed,
                          );                     
    my $not_processed_samples_obj = $all_samples_object->not($samples_processed_obj);
    my @not_processed_samples = @{$not_processed_samples_obj->list};
    print join("\n", @not_processed_samples) . "\n";
}	
SKIP1:
	
############################################################
### Check if event 3 is completed 
### take input of cleaning_threshold
### create collection POOLED_CTX
############################################################
my $exist_pooled_ctx_cols = $db->get_CollectionAdaptor->fetch_by_type("POOLED_CTX");
foreach my $exist_pooled_ctx_col ( @$exist_pooled_ctx_cols ) {
	if ( am_i_the_right_entry($exist_pooled_ctx_col) == 1 ) {
		print "Event 3 completed, POOLED_CTX collection " . $exist_pooled_ctx_col->name . " already exisit, skip create collection\n";
		goto SKIP2;
	}			
}	
### FIXME, how to handle when there are multiple thresholds T
my $complete3 = 0;
my $pooled_ctxs =  $db->get_FileAdaptor->fetch_by_type("POOLED_CTX");
foreach my $pooled_ctx ( @$pooled_ctxs ) {
	next if (am_i_the_right_entry($pooled_ctx) == 0);
	print "Event 3 (pooling ctx) is done\nplease provide cleaning threshold ==> \n";
	my $cleaning_threshold = <>;
	chomp $cleaning_threshold;
	my $collection_name = basename($pooled_ctx->name);
	$collection_name =~ s/\.ctx//;
	$collection_name .= ".user_input_T" . $cleaning_threshold;
	my $collection =  ReseqTrack::Collection->new(
		-name => $collection_name,
		-others => $pooled_ctx,
		-type => "POOLED_CTX",  
		-table_name => 'file',
	);
	throw("Please provide cleanup threshold T") if ( !$cleaning_threshold);
	save_collection($collection, $collection_name, "POOLED_CTX");     		
	$complete3 = 1;
}	
print "Event3 is not completed\n" if ($complete3 == 0);
SKIP2:
			
############################################################			
### Check if event 4 is finished, 
### create collection SAMPLE_CTX_TO_CLEAN for each sample, 
### each collection contains a sample ctx and a cleaned_pool_ctx
############################################################

my $cleaned_pool_ctxs =  $db->get_FileAdaptor->fetch_by_type("POOL_CTX_CLEAN");
if ( !$cleaned_pool_ctxs || @$cleaned_pool_ctxs == 0 ) {
	print "Event 4 is not finished, no file of type POOL_CTX_CLEAN is found\n";
	goto SKIP3;
}

my $exist_sample_ctx_to_clean_cols = $db->get_CollectionAdaptor->fetch_by_type("SAMPLE_CTX_TO_CLEAN");
#print "Number of exist sample_ctx_to_clean collection is " . @$exist_sample_ctx_to_clean_cols . "\n";
my $count = 0;
foreach my $exist_sample_ctx_to_clean_col ( @$exist_sample_ctx_to_clean_cols ) {
	if ( am_i_the_right_entry($exist_sample_ctx_to_clean_col) == 1 ) {
		$count++;
	}
}		
if ($count == $total_sample_cnt) {
	print "Event 4 is finished, collection SAMPLE_CTX_TO_CLEAN is found for all $count samples\n";
	goto SKIP3;
}		

	
my $complete = 0;	
foreach my $cleaned_pool_ctx ( @$cleaned_pool_ctxs ) {
	next if (am_i_the_right_entry($cleaned_pool_ctx) == 0);
	#print "Number of uncleaned sample ctx is " . @uncleaned_sample_ctxs . "\n";
	foreach my $ctx ( @uncleaned_sample_ctxs ) { ### these are per sample uncleaned ctx
		my $collection_name = basename($ctx->name) . ".to_clean";  
		my @array = ($ctx, $cleaned_pool_ctx);
		my $collection =  ReseqTrack::Collection->new(
			-name => $collection_name,
			-others => \@array,
			-type => "SAMPLE_CTX_TO_CLEAN",  
			-table_name => 'file',
		);
		save_collection($collection, $collection_name, "SAMPLE_CTX_TO_CLEAN");			
		$complete = 1;
	}	
}	    						
print "Event 4 is not finished\n" if ($complete == 0);
SKIP3:
				        
############################################################				        
### Check if event 8 (create pooled bubble ctx) is finished
### create collection SAMPLE_CTX_CLEAN_TO_OL_BUBBLE
############################################################
my $bubble_ctxs =  $db->get_FileAdaptor->fetch_by_type("BUBBLE_CTX");
if ( ! $bubble_ctxs || @$bubble_ctxs == 0 ) {
	print "Event 8 is not completed\n";
	goto SKIP4;
}	

my $exist_sample_ctx_clean_to_ol_bubble_cols = $db->get_CollectionAdaptor->fetch_by_type("SAMPLE_CTX_CLEAN_TO_OL_BUBBLE");
my $count1 = 0;
foreach my $sample_ctx_clean_to_ol_bubble_col ( @$exist_sample_ctx_clean_to_ol_bubble_cols ) {
	if (am_i_the_right_entry($sample_ctx_clean_to_ol_bubble_col) == 1) {
		$count1++;
	}
}
if ($count1 == ($total_sample_cnt + 1 )) {
	print "Event 8 is finished, collection SAMPLE_CTX_CLEAN_TO_OL_BUBBLE is found for all $count1 samples (including ref ctx)\n";
	goto SKIP4;
}		
		

my $cleaned_per_sample_ctxs =  $db->get_FileAdaptor->fetch_by_type("CLEANED_PER_SAMPLE_CTX");
my $complete8 = 0;
foreach my $bubble_ctx ( @$bubble_ctxs ) {
	next if  ( am_i_the_right_entry($bubble_ctx) == 0 );
	foreach my $cleaned_per_sample_ctx ( @$cleaned_per_sample_ctxs ) { ### these are per sample cleaned ctx
		next if (am_i_the_right_entry($cleaned_per_sample_ctx) == 0 );
		my $collection_name = basename($cleaned_per_sample_ctx->name) . ".to_bubble";  
		my @array = ($cleaned_per_sample_ctx, $bubble_ctx);
		my $collection =  ReseqTrack::Collection->new(
			-name => $collection_name,
			-others => \@array,
			-type => "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE",  
			-table_name => 'file',
		);
		save_collection($collection, $collection_name, "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE");		
		$complete8 = 1;
	}
	throw ("Please provide ref_binary path") if ( ! $input{ref_binary}); 
	my $ref_binary = $input{ref_binary};
	my $exist_ref_ctx_obj = $db->get_FileAdaptor->fetch_by_name($input{ref_binary});
	my $ref_ctx_objs;
	my @ref_array;
	
	if (!$exist_ref_ctx_obj) {
		$ref_ctx_objs = save_file($ref_binary, "REF_CTX") if ( $input{store} );
		@ref_array = ($ref_ctx_objs->[0], $bubble_ctx);
	}
	else {
		@ref_array = ($exist_ref_ctx_obj, $bubble_ctx);
	}	
	
	my $ref_collection_name = "ref_binary_to_bubble." . $input{pop};
	my $ref_collection =  ReseqTrack::Collection->new(
			-name => $ref_collection_name,
			-others => \@ref_array,
			-type => "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE",  
			-table_name => 'file',
	);
	save_collection($ref_collection, $ref_collection_name, "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE");	
}		
print "Event 8 is not completed\n" if ($complete8 == 0);
SKIP4:


############################################################
### Check if all event 9 (overlay per sample ctx to bubble ctx) finished
### Create collection SAMPLE_CTX_CLEAN_BUBBLED_POOL
############################################################
my %samples_have_bubble_ctx;
my @bubble_ctxs;
my $sample_clean_bubbled_ctxs =  $db->get_FileAdaptor->fetch_by_type("SAMPLE_CTX_CLEAN_BUBBLED");
foreach my $sample_clean_bubbled_ctx ( @$sample_clean_bubbled_ctxs ) {
	next if ( am_i_the_right_entry($sample_clean_bubbled_ctx ) == 0 );
	$samples_have_bubble_ctx{$sample_clean_bubbled_ctx} = 1;
	push @bubble_ctxs, $sample_clean_bubbled_ctx;
}		     


my $exist_sample_ctx_clean_bubbled_pool_cols = $db->get_CollectionAdaptor->fetch_by_type("SAMPLE_CTX_CLEAN_BUBBLED_POOL");
foreach my $sample_ctx_clean_bubbled_pool_col ( @$exist_sample_ctx_clean_bubbled_pool_cols ) {
	if (  am_i_the_right_entry($sample_ctx_clean_bubbled_pool_col) == 1 ) {
		print "Event 9 finished, skip creating collection SAMPLE_CTX_CLEAN_BUBBLED_POOL\n";
		goto SKIP5;
	}
}		
				   
if (keys %samples_have_bubble_ctx == $total_sample_cnt + 1 ) { # ref binary is part of the SAMPLE_CTX_CLEAN_BUBBLED collection
	print "Event 9 (overlay sample ctx with bubble) finished for all samples\n";
	print "**** Please make sure the cortex executable has the correct number of colours!!!*****\n";
	my $collection_name = $input{pop} . "_" . keys(%$samples_to_look_for) . "_bubbled_sample_ctxs_to_pool";
	my $collection =  ReseqTrack::Collection->new(
		-name => $collection_name,
		-others => \@bubble_ctxs,
		-type => "SAMPLE_CTX_CLEAN_BUBBLED_POOL",  
		-table_name => 'file',
	);
	save_collection($collection, $collection_name,"SAMPLE_CTX_CLEAN_BUBBLED_POOL");			
}	
else {
	print "Event 9 not completed for pop $input{pop} yet:\n";
	print "Number of samples to process is " .  keys(%$samples_to_look_for) . "\n";
	print "Number completed so far is " . keys(%samples_have_bubble_ctx) . "\n";
}					        
SKIP5:


############################################################
### Check if event 10 (build multi-colour ctx) finished
### Create collection MULTI_COLOUR_CTX
############################################################
my $multi_colour_ctxs = $db->get_FileAdaptor->fetch_by_type("MULTI_COLOUR_CTX");
if (!$multi_colour_ctxs || @$multi_colour_ctxs == 0 ) {
	print "Event 10 is not finished\n";
	goto SKIP6;
}

my $gts1 = $db->get_FileAdaptor->fetch_by_type("GENOTYPE");
foreach my $gf1 ( @$gts1 ) {
	if (am_i_the_right_entry($gf1) == 1) {
		print "Event 11 finished\n";
		goto SKIP6;
	}
}		

my $bubbles = $db->get_FileAdaptor->fetch_by_type("BUBBLES");
my $correct_bubble;
foreach my $bubble ( @$bubbles ) {
	if ( am_i_the_right_entry($bubble) == 1 ) {
		$correct_bubble = $bubble;
	}
}		
my $complete10 = 0;
foreach my 	$multi_colour_ctx ( @$multi_colour_ctxs ) {
	next if ( am_i_the_right_entry($multi_colour_ctx) == 0 );
	
	my @array = ($multi_colour_ctx, $correct_bubble);
	
	my $collection_name = basename($multi_colour_ctx->name);
	my $collection =  ReseqTrack::Collection->new(
		-name => $collection_name,
		-others => \@array,
		-type => "MULTI_COLOUR_CTX",  
		-table_name => 'file',
	);
	save_collection($collection, $collection_name,"MULTI_COLOUR_CTX");	
	$complete10 = 1;		
}

print "Event 10 is not completed\n" if ($complete10 == 0);				        
SKIP6:	

################################################################################
### Check if event 12 (make covg file from genotype file) is finished
### Carry on creating auxillary table files and create AUX collections (sub-divided)
################################################################################
my $gts = $db->get_FileAdaptor->fetch_by_type("GENOTYPE");
if (!$gts || @$gts == 0 ) {
	print "Event11 (genotyping) is not finished\n";
	goto SKIP9;
}

my $cvgs = $db->get_FileAdaptor->fetch_by_type("COVG_FILE");
if (!$cvgs || @$cvgs == 0 ) {
	print "Event12 (making coverage file from genotype file) has not started\n";
	goto SKIP9;
}
	
my $auxs = $db->get_CollectionAdaptor->fetch_by_type("AUX");
foreach my $aux ( @$auxs ) {
	if ($aux->name =~ /$input{pop}/) {
		print "Event12 (making coverage file from genotype file) is done, auxillary collection " . $aux->name . " with type AUX already exist\n";
		goto SKIP7;
	}
}		

foreach my $gf ( @$gts ) {
	next if (am_i_the_right_entry($gf) == 0);
	my %aux_files_pop;
	my $out_table = create_aux_table_file($input{pop});
	my $covg_file = query_aux_covg_file($gf);
	#print "table file is $out_table\n";
	#print "covg file is $covg_file\n";
	if ( $covg_file eq "null" ) {
		print "Event12 (making coverage file from genotype file) has not finished for " . $gf->name . "\n";
		goto SKIP7;
	}	
		
	if ( -e  $covg_file && -e  $out_table  ) {
		my $table_objs = save_file($out_table, "TABLE_FILE");
		my $covg_objs = save_file($covg_file, "COVG_FILE");		
		
		my @tmp = split(/\_/, basename($covg_file));
		my $sample_ref_number = $tmp[1] + 1;
		
		push @{$aux_files_pop{$input{pop}}}, $covg_objs->[0];
		push @{$aux_files_pop{$input{pop}}}, $table_objs->[0];
						
		my $wc = `wc -l $covg_file`;
		my @tmp3 = split(/\s+/, $wc);
		my $var_cnt = $tmp3[0];
		print "total cnt of var is $var_cnt\n";
		for (my $i=0; $i < $var_cnt/500000; $i++) { 
			my $start = $i * 500000 + 1;					
			my $collection_name = $input{pop} . "_$var_cnt" . "_$sample_ref_number" . "_$start";
			print "collection name is $collection_name\n";
			my $collection =  ReseqTrack::Collection->new(
				-name => $collection_name,
				-others => $aux_files_pop{$input{pop}},
				-type => "AUX",  
				-table_name => 'file',
			);
			save_collection($collection, $collection_name,"AUX") if ($input{run});
		}		
	}
	else {
		print "COVF file $covg_file and/or TABLE file $out_table do not exist\n";
	}	
}	
SKIP7:

#######################################################################
### Check if event 14 (classify results by population) has finished
### Concat the genotype.classified files and store as CLASSIFIED
#######################################################################
my $classified_cols = $db->get_CollectionAdaptor->fetch_by_type("CLASSIFIED_CHUNK");
if ( !$classified_cols || @$classified_cols == 0) {
	print "Event14 classify variant calls hasn't started yet\n";
	goto SKIP8;
}
	
foreach my $classified_col ( @$classified_cols) {
	next if ( $classified_col->name !~ /$input{pop}/);
	my @tmp = split(/\_/,  $classified_col->name);
	my $total_var_num = $tmp[1];
	my $expected_classfied_files = int($total_var_num/500000) + 1;
	print "expected chunky classified files is $expected_classfied_files\n";
	if ( @{$classified_col->others} < $expected_classfied_files) { 
		print "Event14 classify variant calls hasn't finished yet\n";
		print "chunky classified files already made is " . scalar @{$classified_col->others} . "\n";
		goto SKIP8;
	}
	elsif ( @{$classified_col->others} == $expected_classfied_files ) {	
		my $cat_classified_file = $out_dir_pop . "/" . basename($classified_col->name);
		sort_and_cat_files ($classified_col->others, "_", 3, $cat_classified_file);
		if (-e $cat_classified_file ) {
			save_file($cat_classified_file, "CLASSIFIED");	
		}
		else {
			print "catted classified var file does not exist\n";
		}	
	}
	else {
		throw("Number of chunky classified genotype files is more than expected $expected_classfied_files");
	}		
}	

SKIP8:

######################################################################
### Check if event 15 (run_stampy) is finished for all FASTQ files
### SAM_CHUNK collection is completed, if yes
### Create concatecated SAM file and start to generte VCF file
######################################################################
my $fastq_cols = $db->get_CollectionAdaptor->fetch_by_type("5PF_FASTQ");
if (!$fastq_cols || @$fastq_cols ==0 ) {
	print "Event13 (make and split 5pf fastq from genotype file) has not finished yet\n";
	goto SKIP9;
}
		 
my $sam_chunk_cols = $db->get_CollectionAdaptor->fetch_by_type("SAM_CHUNK");
if ( !$sam_chunk_cols || @$sam_chunk_cols ==0 ) {
	print "Event15 run stampy hasn't started\n";
	goto SKIP9;
}
	
foreach my $sam_chunk_col ( @$sam_chunk_cols ) {
	next if ($sam_chunk_col->name !~ /$input{pop}/);
	
	my @tmp = split (/total/, $sam_chunk_col->name);
	my $sam_chunk_cnt = $tmp[1];
	$sam_chunk_cnt =~ s/\.sam//;
	print "Number of chunky sams to cat is $sam_chunk_cnt\n";
	
	if ( @{$sam_chunk_col->others} < $sam_chunk_cnt ) {
		print "Event15 run stampy not finished.\nExpect $sam_chunk_cnt 5pf fastq files to generate sams; have " . scalar(@{$sam_chunk_col->others}) . " sams\n";
		goto SKIP9;
	}
	elsif ( @{$sam_chunk_col->others} == $sam_chunk_cnt	) {		
		
		my $genotype_files = $db->get_FileAdaptor->fetch_by_type("GENOTYPE");
		my $geno_basename;
		my $num = 0;
		foreach my $genotype_file ( @$genotype_files ) {
			if ($genotype_file->name =~ $input{pop}) {
				$geno_basename = basename($genotype_file->name);
				throw("Genotype file " . $genotype_file->name . " does not exist") unless ( -e $genotype_file->name);
				$num++;
			}
		}	
		throw("No genotype file found") if ($num ==0);
		my $cat_sam_file = $out_dir_pop . "/" . $geno_basename . ".5pflanks.sam";
		throw("catted sam $cat_sam_file already exist") if (-e $cat_sam_file);
		sort_and_cat_files ($sam_chunk_col->others, ".", 1, $cat_sam_file);
		if (-e $cat_sam_file ) {
			my $cat_sam_objs = save_file($cat_sam_file, "SAM");
			my $collection_name = $sam_chunk_col->name;
			print "collection name is $collection_name\n";
			my $collection =  ReseqTrack::Collection->new(
				-name => $collection_name,
				-others => $cat_sam_objs,
				-type => "SAM",  
				-table_name => 'file',
			);
			save_collection($collection, $collection_name,"SAM");		
		}
		else {
			print "Catted sam file 	$cat_sam_file does not exist\n";
		}		
	}
	else {
		throw("Number of chunky sam files is more than expected $sam_chunk_cnt");
	}							
}


SKIP9:


print "End\n";

##############
#### SUBS ####
##############
sub sort_and_cat_files {
	my ($others, $div, $pos_in_file, $out_cat_file) = @_;
	
	warning("Concatenated file $out_cat_file already exists") if (-e $out_cat_file);
	
	my %order_files;	
	foreach my $file ( @$others ) {	
		throw("File to cat " . $file->name . "does not exist") unless ( -e $file->name);	
		my @tmp;
		if ($div eq "." ) {
			@tmp = split(/\./, basename($file->name) );
		}
		elsif ( $div eq "_" ) {	
			@tmp = split(/\_/, basename($file->name) );
		}
		my $index = $tmp[$pos_in_file];	
		$order_files{$index} = $file->name;
	}	
			
	my @sorted_index = sort {$a<=>$b} ( keys %order_files );
	`grep ^@ $order_files{$sorted_index[0]} > $out_cat_file` if $input{run}; 
	## to clean up output cat file as well.  
	## For sam file to get header, classified files do not have header
		
	foreach my $j ( @sorted_index ) {
		print "Cat number $j file: $order_files{$j} into $out_cat_file.... if -run\n";
		`cat $order_files{$j} | grep -v ^@ >> $out_cat_file` if $input{run};
		my $exit = $?>>8;
		throw("cat files failed for $order_files{$j}\n") if ($exit >=1);	
	}	
	print "Catted file is $out_cat_file\n";
	
	return 1;
}

		
sub create_aux_table_file {
	my ($pop) = @_;
	
	my $genotype_logs = `find $input{log_dir} -name \"\*$pop*genotype\*out\" -print`;
	my @logs = split(/\n/, $genotype_logs);
	my $cnt = 0;
	my $genotype_log;
	foreach my $log ( @logs ) {
		next if $log !~ /event11/;
		#print $log . "\n";
		$cnt++;
		$genotype_log = $log;
	}	

	throw("More than one genotype log files found for $input{pop} in $input{log_dir}") if $cnt > 1;
	
	my $log_base = basename($genotype_log);
	
	my $output_table = $out_dir_pop  . "/" . $log_base . ".table";
	
	my $command = "perl /nfs/1000g-work/G1K/work/bin/cortex/scripts/analyse_variants/make_read_len_and_total_seq_table.pl $genotype_log >& $output_table";
	print "Running ...\n$command if -run\n";
	
	`$command` if ($input{run});
	my $exit = $?>>8;
	throw("make_read_len_and_total_seq_table.pl failed for $genotype_log\n") if ($exit >=1);	
	return 	$output_table;
}


sub query_aux_covg_file {
	my ($gf) = @_;
	my $covg_file_name = $gf->name . ".covg_for_classifier";
	my $cvg = $db->get_FileAdaptor->fetch_by_name($covg_file_name);
	if ( $cvg ) {
		return $cvg->name;
	}
	else {
		return "null";
	}	
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


sub am_i_the_right_entry {  ### FIXME, need to check for other parts in the file name to make sure that they match as well.
	my ($entry) = @_;
	my $flag = 0;
	#print "entry is " . $entry->name . "\n";
	my @tmp = split(/\./, basename($entry->name));
	if ( 	$entry->type eq "SAMPLE_CTX" || 
			$entry->type eq "CLEANED_PER_SAMPLE_CTX" || 
			$entry->type eq "SAMPLE_CTX_TO_CLEAN" ||
			$entry->type eq "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE" ||
			$entry->type eq "SAMPLE_CTX_CLEAN_BUBBLED" ) {
		
		if ( $entry->name =~ /ref|human/i ) { ### FIXME, need better ways to identify reference 
			if ( $entry->type eq "SAMPLE_CTX_CLEAN_TO_OL_BUBBLE" || $entry->type eq "SAMPLE_CTX_CLEAN_BUBBLED") {
				if ($entry->name =~ /$input{pop}/ ) {
					$flag = 1;
				}
			}
			else {
				$flag = 1;
			}
		}	
		elsif ( $sample_to_pop_mapping->{$tmp[0]} eq $input{pop} ) { ### FIXME, need to check for other parts in the file name to make sure that they match as well.
			$flag = 1;
		}
	}
	else {	
		if ( $entry->name =~ /$input{pop}/ ) { ### FIXME, need to check for other parts in the file name to make sure that they match as well.
			my @tmp = split (/\_/, 	basename($entry->name) );
			my $sample_cnt = $tmp[1];
			if ($sample_cnt eq "cleaned") { ### This is to handle the pool_ctx_clean that have been renamed and put up to the ftp site (rerun for GIH, YRI, LWK)
				$flag = 1;
			}	
			elsif ($sample_cnt == $total_sample_cnt ) {
				$flag = 1;
			}
		}
	}	    	        
	return $flag;
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
=head1 NAME

ReseqTrack/scripts/process/cortex_progress_manager.pl

=head1 SYNOPSIS

This script is the driver of run cortex pipeline.  It runs periodically to check the status of each job and the existence of appropriate output
files from each job. Once all jobs for a given event are accomplished, the progress manager would create appropriate collection in the db that will 
serve as input to kick off downstream events.

After event 11 (genotype) is done, the pipeline enters the "post processing" stage. Run the progress manager with a -post_processing tag.

=head1 EXAMPLES    

Regular processing:
perl $ZHENG_RB_VC/scripts/process/cortex_process_manager.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index  \
-ref_binary /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37.proper_chroms.k31.ctx \
-pop LWK \
-run -store -save_collection -update \

Post processing:
perl $ZHENG_RB_VC/scripts/process/cortex_process_manager.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-seq_index /nfs/1000g-work/G1K/work/zheng/cortex/20120522.sequence.index.ori \
-ref_binary /nfs/1000g-archive/vol1/ftp/technical/working/20120814_cortex_resources/human_g1k_v37.proper_chroms.k31.ctx \
-post_processing \
-log_dir /nfs/nobackup/resequencing_informatics/zheng/run_cortex/log \
-out_dir /nfs/1000g-work/G1K/work/zheng/cortex/post_process \
-pop PEL \
-run -store -save_collection -update \
