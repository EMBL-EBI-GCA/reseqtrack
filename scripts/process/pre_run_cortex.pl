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
	%sample_se,
	%sample_pe1,
	%sample_pe2,
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
	'sample_list:s',
	'sample:s',
		
  	'output_dir=s',
  	'output_file_type:s',

  	'store!',
  	'save_collection!',
  	'run!',
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
$input{output_file_type} = 'FASTQ_LIST' if (!$input{output_file_type} );
   
if (!$input{seq_index}  && !$input{pop}) {
	throw("Please provide a sequence index file and a population name to be processed");
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
);

my $fa = $db->get_FileAdaptor;

my $lines = get_lines_from_file($input{seq_index});
throw("Seq index file $input{seq_index} does not exist") if (!$lines || @$lines==0);

my %samples_to_process;
if (defined $input{sample_list} ) {
	my $samples = get_lines_from_file($input{sample_list}) ;
	foreach my $sample (@$samples) {
		$samples_to_process{$sample} = 1;
	}	
}
elsif ( defined $input{sample} ) {
	$samples_to_process{$input{sample}} = 1;
}	
	
foreach my $line (@$lines) {
	next if ( $line =~ /FASTQ_FILE/i ); ### skip the headline
	my @bits = split(/\t/, $line);
	my $sample_name = $bits[9];	
	my $fastq = $bits[0];
	my $platform = $bits[12];
	my $analysis_grp = $bits[25];
	my $withdrawn = $bits[20];
	my $pop = $bits[10];
		
	if ($withdrawn != 1 && 
		$platform =~ /ILLUMINA/i && 
		$analysis_grp =~ /low coverage/) {
			next if ( defined $input{pop} && $pop ne $input{pop} );
			if ($input{sample_list} || $input{sample}) {
				if ($samples_to_process{$sample_name}) {
					print "input fastq is $fastq\n";
					process_cortext_input($fastq, $sample_name); ### populate a hash
				}
			}	
			else {
				print "input fastq is $fastq\n";
				process_cortext_input($fastq, $sample_name);
			}	
	} 
}

print_and_load_list(\%sample_se, "se");
print_and_load_list(\%sample_pe1, "pe1");
print_and_load_list(\%sample_pe2, "pe2");

######### SUBS #########
sub print_and_load_list {
	my ($hash, $type) = @_;
	foreach my $sam ( keys %$hash ) {	
		my $output_list = $input{output_dir} . "/" . $sam . ".$type" . "_list";		
		open (OUT, ">", $output_list) || throw("Cannot open output file $output_list");
		print OUT join("\n", @{$hash->{$sam}} ) . "\n";
		load_fastq_list($output_list, $sam);
	}
	return 1;
}

sub load_fastq_list {
	my ($file, $sample) = @_;
	
	my $loader = ReseqTrack::Tools::Loader::File->new
		(
		   -file => [$file],
		   -do_md5 => 1,
		   -size => (-s $file),
		   -hostname => $input{host},
		   -db => $db,
		   -assign_types => 0,
		   -check_types => 0,
		   -type => $input{output_file_type},
		   -update_existing => $input{update},
		);
		
	my $file_objects;
	if($input{store}){
		  $loader->process_input();
		  $loader->create_objects();
		  $loader->sanity_check_objects();
		  $file_objects = $loader->load_objects();
	}	
	
	if ( $input{save_collection} ) {
		my $collection_name = $sample;
		my $collection =  ReseqTrack::Collection->new(
			-name => $collection_name,
			-others => $file_objects,
			-type => $input{output_file_type},  
			-table_name => 'file',
		);
		$db->get_CollectionAdaptor->store($collection);  ## new files will be added in without over write the old ones.
	}
		
	return 1;
}

sub process_cortext_input {
	my ($fq, $samp) = @_;
	
	my $basename = basename($fq);
	
	my $outdir = "$input{output_dir}/$samp";
	mkpath($outdir) unless (-e $outdir);
	my $old_path = "/nfs/1000g-archive/vol1/ftp/" . $fq;
	my $new_path = $outdir . "/" . $basename;
	print "rsync $old_path $outdir\n";
	eval {
		`rsync $old_path $outdir/` if ($input{run});
	};
	throw("rsync file $old_path failed, $@\n") if $@;

=head	
	eval {
		`gunzip $new_path` if ($input{run});
	};	
	throw("gunzip file $new_path failed, $@\n") if $@;	
	
	$new_path =~ s/.gz//;
=cut		
	if ($basename =~ /_1.filt.fastq/) {
		push @{$sample_pe1{$samp}}, $new_path;
	}
	elsif (	$basename =~ /_2.filt.fastq/) {
		push @{$sample_pe2{$samp}}, $new_path;
	}
	elsif ( $basename =~ /filt.fastq/) {
		push @{$sample_se{$samp}}, $new_path;
	}
	else {
		throw("File $basename is not a fastq file");
	}
	
	return 1;
}		 
	
	
=pod


perl $ZHENG_RT/bin/pre_run_cortex.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index \
-pop LWK \
-sample_list /nfs/1000g-work/G1K/work/zheng/cortex/LWK_test/sample_list	\
-store \
-save_collection \
-update \


__DATA__

1000genomes.ebi.ac.uk> perl $ZHENG_RB_VC/scripts/process/pre_run_cortex.pl -seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index -pop LWK |wc -l
1115
1000genomes.ebi.ac.uk> perl $ZHENG_RP/bin/run_cortex.pl -seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index -pop ASW | wc -l
844
1000genomes.ebi.ac.uk> perl $ZHENG_RP/bin/run_cortex.pl -seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index -pop GIH | wc -l
3525
1000genomes.ebi.ac.uk> perl $ZHENG_RP/bin/run_cortex.pl -seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index -pop YRI | wc -l
4051





