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
use ReseqTrack::Tools::ConvertBam;
use ReseqTrack::Tools::myTIME;
use ReseqTrack::Tools::Loader::Archive;

#my $start_time = time();

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

	'bam=s',
	'collection=s',
	'collection_type=s',
	
	'program=s',
	'parameters=s',
	'reference_fasta=s',

  	'output_dir=s',
  	'output_file_type:s',

  	'store!',
  	'archive!',
 	'update!',
 	'host:s',
  	'help!',		
);

$input{update} = 0 if (!$input{update});

if (!$input{output_dir} ) {
	$input{output_dir} = `pwd`;
	chomp $input{output_dir};
}

$input{output_dir} =~ s/\/$//;
	
$input{host} = '1000genomes.ebi.ac.uk' if (!$input{host});
$input{program} = '/nfs/1000g-work/G1K/work/bin/crammer' if (!$input{program});
$input{output_file_type} = "CRAM" if ( !$input{output_file_type});

if (!$input{output_file_type} && $input{store} ) {
	throw("Please provide a file type if you want to store the output file in the database");
}	

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
);

my $bam;

if ( 	$input{collection} && $input{collection_type} ) {
	my $collection =
	    $db->get_CollectionAdaptor->fetch_by_name_and_type($input{collection},
	                                                        $input{collection_type});
	throw("No collection found for $input{collection} and $input{collection_type}") if (!$collection); 
	
	foreach my $list (@{$collection->others}){
	    print "list: " . $list->name . "\n";
	    if ($list->name =~ /chrom20/) {
	        $bam = $list->name;
	    }    
	}
}
elsif ( $input{bam} ) {
	$bam = $input{bam};
}	
else {
	throw("Please input file, either as collection name and type, or as bma path");
}	

throw("Bam file $bam does not exist") unless (-e $bam);

$db->dbc->disconnect_when_inactive(1); 

my $parameter_hash = parameters_hash( $input{parameters} );

my $output_file_name = assign_cram_name($bam);  

my $cram = ReseqTrack::Tools::ConvertBam->new (
	-input_files => $bam,
	-output_format => 'cram',
	-options => $parameter_hash,
	-program => $input{program},
	-reference_fasta => $input{reference_fasta},
	-output_file => $output_file_name,
);

$cram->run;

my ($output_cram) = @{$cram->output_files};

print "output file is $output_cram\n";
my $crai = $output_cram . ".crai";
save_file($output_cram, "CRAM") if $input{store};

if ( -e $crai) {
	save_file($crai, "CRAI") if $input{store};
}
else {
	throw("Crai file $crai does not exist");
}	

archive_files($output_cram, $crai) if ($input{archive} &&  -e $crai );

##### SUBS ####
sub archive_files {
	my ($f1, $f2) = @_;
	
	my ($time_stamp, $month_stamp, $day_stamp) = ReseqTrack::Tools::myTIME::get_time();
	my $cram_basename = basename($f1);
	
	my $archive_list = '/nfs/1000g-work/G1K/scratch/zheng/tmp/' . $cram_basename . ".tmp_archive_list." . $time_stamp;

	open (LIST, ">", $archive_list) || throw("Cannot open temparary archive list $archive_list\n");
	print LIST "$f1\n$f2\n";
	close LIST;
	
	my $action_string = "archive";

	my $max_number = 1000;
	my $priority = 50;
	my $verbose = 0;
	
	my $archiver = ReseqTrack::Tools::Loader::Archive->new(
                                                       -list_file => $archive_list,
                                                       -dbhost => $input{dbhost},
                                                       -dbname => $input{dbname},
                                                       -dbuser  => $input{dbuser},
                                                       -dbpass  => $input{dbpass},
                                                       -dbport  => $input{dbport},
                                                       -action => $action_string,
                                                       -verbose => $verbose,
                                                       -priority=> $priority,
                                                       -max_number=>$max_number,
                                                       -no_lock => 1,
                                                      );


	$archiver->process_input();
	#$archiver->cleanup_archive_table($verbose); #leave this out equals to have the -skip_cleanup option
	$archiver->sanity_check_objects();
	$archiver->archive_objects() if $input{archive};
	return 1;
}

	
sub assign_cram_name {
	my ($bam) = @_;
	my $output;
	my $bam_basename =  basename($bam);
	my ($sample) = split(/\./, $bam_basename);
	my $dir;
	#print "bam basename is $bam_basename\n";
	if ($bam_basename =~ /low_coverage/) {
        $dir = $input{output_dir} . "/$sample/alignment/";
        mkpath($dir) unless (-e $dir);
        $output = $dir . $bam_basename . ".cram";
    }
    elsif ( $bam_basename =~ /exome/) {   
        $dir = $input{output_dir} . "/$sample/exome_alignment/";
        mkpath($dir) unless (-e $dir);
        $output = $dir . $bam_basename . ".cram";
    }
    else {
		throw("bam file $bam_basename is neither low coverage nor exome");
    }	
    #print "output file is $output\n";
    return $output;
}

    
sub parameters_hash{
  my ($string) = @_;

  my %parameters_hash;
  print "parameter string is $string\n";
	
  if ($string) {
    
    my @pairs;
    if ($string !~ /,/  ) {
                push @pairs, $string;
    }
    elsif($string =~  /,/ && ($string =~ /:/ || $string =~ /=>/) ){
                @pairs = split (/,/, $string);
    }
    else {
          throw("Please provide running parameters as name and value separated by ':' or '=>';  Multiple name and value pairs can be used and should be separated by ',' ");
    }

    foreach my $pair(@pairs){
                my ($key, $value);
                if ($pair =~ /:/) {
                   ($key, $value) = split (/:/, $pair); 
                }
                elsif ( $pair =~ /=>/) {
                    ($key, $value) = split (/=>/, $pair); 
                }
                else {
                    throw("Please provide running parameters as name and value separated by ':' or '=>'");
                }    
                $key   =~ s/^\s+|\s+$//g;
                $value =~ s/^\s+|\s+$//g;
                $parameters_hash{$key} = $value;
                print "key is $key, value is $value\n";
    }#end of foreach pairs
  }# end of if ($string)
 
  return \%parameters_hash;
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
perl /nfs/1000g-work/G1K/work/zheng/reseqtrack_branch_zheng/crammer/scripts/process/run_crammer.pl -dbhost mysql-g1kdcc-public -dbuser g1krw -dbpass thousandgenomes -dbport 4197 -dbname g1k_archive_staging_track \
-program /nfs/1000g-work/G1K/work/bin/crammer \
-collection HG01377.ILLUMINA.bwa.low_coverage \
-collection_type BAM \
-output_dir /tmp \
-reference_fasta /nfs/1000g-archive/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz \
-parameters "capture-all-tags=>1,preserve-read-names=>1,heap_size=>4g,lossy-quality-score-spec=>m999,create-index=>1" \
-store \
-output_file_type CRAM

OR

-bam /nfs/1000g-archive/vol1/ftp/data/HG01377/alignment/HG01377.chrom20.ILLUMINA.bwa.CLM.low_coverage.20120522.bam \