#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use Getopt::Long;

my (
    $dbhost,
    $dbuser,
    $dbpass,
    $dbport,
    $dbname,
    $file,
    $chunk_size,
    $chunk_type,
    $run,
    );
    
    
&GetOptions( 
  'dbhost=s'      		=> \$dbhost,
  'dbname=s'      		=> \$dbname,
  'dbuser=s'      		=> \$dbuser,
  'dbpass=s'      		=> \$dbpass,
  'dbport=s'      		=> \$dbport,
	
  'file=s'				=> \$file,  #### input file can be a BAM header that contains the reference assembly information
  'chunk_size=i'		=> \$chunk_size,
  'chunk_type=s'		=> \$chunk_type,
  'run!'				=> \$run,
  );
  

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );
    
my $isa = $db->get_InputStringAdaptor();
       
open (IN, "<", $file) || throw("Cannot open chrom length file $file");

while (<IN>) {
	chomp;
	my @tmp = split(/\t/, $_);
	my ($foo, $chrom) = split(/:/, $tmp[1]);
	my ($foo2, $length) = split(/:/, $tmp[2]);
	
	for (my $i = 0; $i < $length/$chunk_size; $i++) {
		
		my $chunk_name;
		
		if ( ($i+1)*$chunk_size < $length) {
			$chunk_name = $chrom . ":" . ($i*$chunk_size +1) . "-" . ($i+1)*$chunk_size;
		}
		else {
			$chunk_name = $chrom . ":" . ($i*$chunk_size +1) . "-" . $length;
		}
					
		my $is_obj = ReseqTrack::InputString->new(
               -name => $chunk_name,
               -type => $chunk_type,
             );
       
       print "name: $chunk_name\n";
             
       $isa->store($is_obj) if ($run);
	}
}	             

=pod	
This script is to parse a BAM header file, get the chromosome length for each chromosome and store in the input string table chunks of user-defined length.

Example:
perl insert_chrom_chunk_input_string_tb.pl $WRITE_DB_ARGS -dbname zheng_var_call -file /nfs/1000g-work/G1K/work/zheng/snp_calling/chrom_length_ncbi37 -chunk_size 3000000 -chunk_type CHR_CHUNK	
