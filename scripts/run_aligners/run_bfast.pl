#!/usr/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::RunAlignment::BFAST;
use Getopt::Long;
use Data::Dumper;

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4175;
my $dbname;
my $run_id;
my $type;
my $output_dir;
my $build;
my $read_length;
my $working_dir;
&GetOptions( 
	    'dbhost=s'      => \$dbhost,
	    'dbname=s'      => \$dbname,
	    'dbuser=s'      => \$dbuser,
	    'dbpass=s'      => \$dbpass,
	    'dbport=s'      => \$dbport,
	    'run_id=s' => \$run_id,
	    'type=s' => \$type,
	    'output_dir=s' => \$output_dir,
	    'build=s' => \$build,
	    'read_length=s'=>\$read_length,
	   );



my $db = ReseqTrack::DBSQL::DBAdaptor->new(
					   -host   => $dbhost,
					   -user   => $dbuser,
					   -port   => $dbport,
					   -dbname => $dbname,
					   -pass   => $dbpass,
					  );

my $ca = $db->get_CollectionAdaptor;

my $collection = $ca->fetch_by_name_and_type($run_id, $type);

throw("Failed to find a collection for ".$run_id."      ".$type." from ".  $dbname) 
  unless($collection);

my $run_alignment = ReseqTrack::Tools::RunAlignment::BFAST->new(
								-input => $collection,
								-working_dir => $output_dir,
								-name => $run_id,
								-read_length=>$read_length,
								-working_dir=>'/nfs/1000g-work/G1K/scratch/smithre/BFAST_TEST',
							       );
    
$run_alignment->create_cmds;
$run_alignment->run();
  
   

