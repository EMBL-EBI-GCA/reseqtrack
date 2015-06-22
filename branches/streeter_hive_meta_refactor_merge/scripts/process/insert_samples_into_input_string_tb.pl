#!/usr/bin/env perl

use strict;
use warnings;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use Getopt::Long;
use ReseqTrack::Tools::FileSystemUtils;
use File::Basename;

my %input;
    
&GetOptions( 
  \%input,
  'dbhost=s',
  'dbname=s',
  'dbuser=s',
  'dbpass=s',
  'dbport=s',	
  'seq_index=s',
  'pop=s',
  'run!',
  );
  

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $input{dbhost},
  -user   => $input{dbuser},
  -port   => $input{dbport},
  -dbname => $input{dbname},
  -pass   => $input{dbpass},
    );
    
my $isa = $db->get_InputStringAdaptor();

my $lines = get_lines_from_file($input{seq_index});
throw("Seq index file $input{seq_index} does not exist") if (!$lines || @$lines==0);

my @tmp = split(/\./, basename($input{seq_index}));
my $date;
foreach my $piece ( @tmp ) {
	if ( $piece =~ /\d+/) {
		$date = $piece;
	}
}

my $input_string_type = $input{pop} . "_" . $date;	
my %samples_to_insert;	
foreach my $line ( @$lines ) {
	my @bits = split(/\t/, $line);
	my $sample_name = $bits[9];	
	my $fastq = $bits[0];
	my $platform = $bits[12];
	my $analysis_grp = $bits[25];
	my $withdrawn = $bits[20];
	my $pop = $bits[10];
		
	if ($pop eq $input{pop} &&
		$withdrawn != 1 && 
		$platform =~ /ILLUMINA/i && 
		$analysis_grp =~ /low coverage/) {
		
		$samples_to_insert{$sample_name} = 1;
	}	
}	

foreach my $s ( keys %samples_to_insert ) {
							
		my $is_obj = ReseqTrack::InputString->new(
               -name => $s,
               -type => $input_string_type,
             );
       
       print "name: $s\t$input_string_type\n";
             
       $isa->store($is_obj) if ($input{run});
}	             

=pod	
perl $ZHENG_RB_VC/scripts/process/insert_samples_into_input_string_tb.pl $WRITE_DB_ARGS -dbname zheng_run_cortex \
-seq_index /nfs/1000g-archive/vol1/ftp/sequence_indices/20120522.sequence.index \
-pop LWK \
-run
