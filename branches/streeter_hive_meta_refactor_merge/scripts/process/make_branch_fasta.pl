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
	'seq_index=s',
	
	'kmer_size=i',
	
  	'output_dir=s',
  	'output_file_type:s',

  	'store!',
  	'save_collection!',
 	'update!',
 	'host:s',
  	'help!',		
);


$input{update} = 0 if (!$input{update});

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd` ;
	chomp $input{output_dir};
}

$input{host} = '1000genomes.ebi.ac.uk' if (!$input{host});
   
my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
);

my @inputs;
if ( 	$input{collection} && $input{collection_type} ) {
	my $collection =
	    $db->get_CollectionAdaptor->fetch_by_name_and_type($input{collection},
	                                                        $input{collection_type});
	throw("No collection found for $input{collection} and $input{collection_type}") if (!$collection); 
	
	foreach my $list (@{$collection->others}){
	    print "list: " . $list->name . "\n";
	    push @inputs, $list->name;
	}
	my ($samples_to_look_for, $sample_to_pop_mapping) = parse_seq_index($input{seq_index});
	$pop = $sample_to_pop_mapping->{$input{collection}};
}

$db->dbc->disconnect_when_inactive(1); 

if ( @inputs > 1) {
	throw ("There are two bubbles associated with collection $input{collection} and $input{collection_type}");
}

my $bubble_file = $inputs[0];

my $outfile = $input{output_dir} . "/$pop/" . basename($bubble_file) . ".branches.fasta";

`rm $outfile` if (-e $outfile);

my $cmd = "/nfs/1000g-work/G1K/work/bin/cortex/scripts/analyse_variants/make_branch_fasta.pl --callfile $bubble_file --kmer " . $input{kmer_size};
print "command line is \n$cmd\n";
system($cmd);

my $exit = $?>>8;
throw("make_branch_fasta failed for collection $input{collection}\n") if ($exit >=1);

$db->dbc->disconnect_when_inactive(1); 

if (-e $outfile ) {
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


##############
#### SUBS ####
##############
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
perl $ZHENG_RB_VC/scripts/process/make_branch_fasta.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-collection LWK_2_samples.k31.t1.bubbles \
-collection_type BUBBLES \
-store \
-output_file_type BUBBLE_FASTA \
-save_collection \
-kmer_size 31	\
-seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index	
