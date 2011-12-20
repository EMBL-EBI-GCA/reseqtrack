#!/usr/local/bin/perl
use warnings;
use strict;

use Getopt::Long;
use ReseqTrack::Tools::Exception;
use ReseqTrack::Tools::FileUtils;
use ReseqTrack::Tools::FileSystemUtils;
use ReseqTrack::Tools::SequenceIndexUtils;
use ReseqTrack::Host;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::DBSQL::FileAdaptor;
use File::Basename;

use Data::Dumper;

my %input;
my %snps;
my $snps_file = "/nfs/1000g-work/G1K/work/REFERENCE/snps/grc37_snps/OMNI_1525_PASS_V2/bin/20110601_omni1525_v2.bin.list";

#eclipse test. smithre

$input{dbname} = "g1k_archive_staging_track";
$input{dbport} = 4197;
$input{dbhost} = "mysql-g1kdcc-public";
$input{dbuser} = "g1kro";

GetOptions(
	\%input,    'dbhost=s', 'dbname=s',      'dbuser=s',
	'dbpass=s', 'dbport=s', 'working_dir=s', 'verbose!',
	   'skip_name_check!',
);





my $db = &get_db_connection( \%input );



print "\n\nStarting Jr Health check of: " . $input{dbname} . "\n";
print '-' x 45;
print "\n";

open my $SNPS ,'<', $snps_file || die "Nope: snps_file\n";
while (<$SNPS>){
  chomp;
  $snps{$_} = 1;
}
close $SNPS;
&collections_with_missing_file_objects($db);

&failed_glf_but_FILTERED_FASTQ ($db);

&history_objects_with_no_assoc_obj( $db, 'file');

&history_objects_with_no_assoc_obj( $db, 'collection');

&objects_with_no_history ( $db, 'collection');

&fastq_attrib_stat_obj_check( $db, 'FILTERED_FASTQ', $input{verbose} );

&fastq_attrib_stat_obj_check( $db,'WITHDRAWN_FILTERED_FASTQ',$input{verbose});


&sample_swap_results ( $db, $input{verbose});

&more_read_than_bases( $db, $input{verbose} );

&no_genotype_results_for_run_id( $db, $input{verbose} );


exit if $input{skip_name_check};

&check_fastq_collection_types_consistent( $db, 'FILTERED_FASTQ',
	$input{verbose} );

##&check_fastq_collection_types_consistent( $db, 'WITHDRAWN_FILTERED_FASTQ',
##	$input{verbose} );

print "\n\nDone \n";

sub collections_with_missing_file_objects {
 my ( $db,$verbose) = shift;

my $sql ="select cg.*, c.name, c.type from collection_group cg,collection c   where cg.other_id not in (select file_id from file) and cg.collection_id not in (55697,55698,55699,55700) and cg.collection_id = c.collection_id" ;

 print "Looking for collections with missing file_ids linked to them\n\n";

 print $sql,"\n\n" if $verbose;
   
 my $sth = $db->dbc->prepare($sql);
 $sth->execute;
 my $ctr = 0;
 
 print "coll_id cg_other_id coll_name  coll_type\n";
 print "======= ==========  =========  =========\n";
 while ( my @row = $sth->fetchrow_array ) {
   $ctr++;
   my $line = join ( "\t", @row);
   print $line ,"\n";
 }

 print "Found $ctr collections with missing file_id linked to them\n";

  $sth->finish;
 return;

}





=head2 history_objects_with_no_coll_obj 

  Arg [1]   : DB adaptor
  Arg [2]   : verbose
 
  Function  : collections get deleted, but not there history obj

=cut

sub history_objects_with_no_assoc_obj {
    my ( $db, $table,$verbose ) = @_;
    
    my $get = "select other_id,table_name  from history where table_name = \"$table\" ";
    my $clause ="  other_id not in (select ${table}_id from $table) group by other_id,table_name"; 

    my $sql = "$get AND  $clause" ;
    
    print $sql,"\n\n" if $verbose;
   
    my $sth = $db->dbc->prepare($sql);
    $sth->execute;
    my $ctr = 0;

    while ( my @row = $sth->fetchrow_array ) {
            $ctr++;
      }
    $sth->finish;

 print "Got $ctr  history (table = $table) objects with no existing $table objects\n";
    return;
}


=head2




=cut

sub failed_glf_but_FILTERED_FASTQ {
	my ( $db, $verbose ) = @_;
	my  @results;

	print "Checking for genotype FAILED  results\n";
	
	my $sql = "select name from genotype_results where verdict = \"FAILED\"";
	my $sth = $db->dbc->prepare($sql);
	$sth->execute;


	while ( my $rowArrayref = $sth->fetchrow_arrayref ) {
        
		#print "$ctr @$rowArrayref[0] @$rowArrayref[1]\n";
		push( @results, @$rowArrayref[0] );
	}
	print "Got ", scalar(@results), " runs that FAILED genotype check\n";

	my $ca = $db->get_CollectionAdaptor;

	foreach my $name (@results){
	  my $coll = $ca->fetch_by_name_and_type ($name,"WITHDRAWN_FILTERED_FASTQ");

	 if (! $coll){
	   print "Possible error: Got $name failed GLF check but has no ";
	   print "WITHDRAWN_FILTERED_FASTQ collection\n";

	  }

	} 

 return;
}

=head2 objects_with_no_history 

  Arg [1]   : DB adaptor
  Arg [2]   : verbose
 
  Function  : Any files with no History objects


=cut

sub objects_with_no_history {
    my ( $db, $table,$verbose ) = @_;
    my $scratch ='scratch/staging_area/sequence_staging';
    #print "Any $table objects with no History objects\n";
    
    my $sql =
      "select ${table}_id, name from ${table} where ${table}_id not in "
      .  "(select other_id from history where table_name = \"$table\") ";
    
    print $sql,"\n\n" if $verbose;
    
    my $sth = $db->dbc->prepare($sql);
    $sth->execute;
    my $ctr = 0;

    while ( my @row = $sth->fetchrow_array ) {
    	if ($row[1] =~ /scratch/){
    	#	next;
    	}
        if ($row[1] =~ /withdrawn/){
         #   next;
        }
    	 if ($row[1] =~ /pilot/){
          #  next;
        }
         if ($row[1] =~ /sequence_read/){
           # next;
        }
        
        print "@row\n" if $verbose;
        $ctr++;
    }
    
    print "Have $ctr case(s) of $table objects where no history objects present\n";
    
    $sth->finish;



    return;
}





=head2 sample_swap_results

  Arg [1]   : DB adaptor
  Arg [2]   : verbose
 
  Function  : See if missed re-running any sample swaps


=cut

sub sample_swap_results {
	my ( $db, $verbose ) = @_;
	
	#print "Checking for genotype sample swap results\n";
	
	my $sql =
	  'select gt.name, gt.claimed , rmi.run_id, rmi.sample_name from genotype_results gt , run_meta_info rmi ' .
	    'where gt.name = rmi.run_id and gt.claimed != rmi.sample_name';
  
   # print $sql,"\n";
	
	my $sth = $db->dbc->prepare($sql);
	$sth->execute;
	my $ctr = 0;

	while ( my @row = $sth->fetchrow_array ) {
		print "@row\n" if $verbose;
		$ctr++;
	}
	
	print "Have $ctr case(s) where genotype_results.claimed != rmi.sample_name\n";
	
	$sth->finish;



	return;
}



=head2 check_withdrawn_collection_types

  Arg [1]   : DB adaptor
  Arg [2]   : verbose
 
  Function  : Occasionally stuff is withdrawn and type are not changed.


=cut

sub check_fastq_collection_types_consistent {
	my ( $db, $type, $verbose ) = @_;
	my $check;

	print "Checking $type collection file types\n";

	if ( $type =~ /WITHDRAWN/ ) {
		
		$check = '/vol1/ftp';
	}
	else {
		$check = 'withdrawn';
	}
	
    print "checking for $check in $type file names\n";
	
	my $ca = $db->get_CollectionAdaptor;

	my $collections = $ca->fetch_by_type("$type");

	throw("No collection object is found; please check the run_id and type\n")
	  if ( !$collections );

	print "Have ", scalar(@$collections), " collections\n";
    my $ctr = 0;
    my %tracker;
	foreach my $coll (@$collections) {
		my $others = $coll->others;
		foreach my $q (@$others) {
			if (  $q->name =~ /withdrawn/g  ) {
				print $coll->name, "\t", $coll->type, "\t", $q->name, "\n";
				push (@{$tracker{$coll->name}}, $q->name);
			}
		}
	}
	
	my $bad_collections = scalar ( keys %tracker);
	print "Got $bad_collections bad collections with mixed types\n";
	
#	foreach my $key (keys %tracker){
#		my @arr = $tracker{$key};
#		print $key,"\n";
#		foreach (@arr){
#		  print;
#		  print "\n";	
#		}		
#	}
	
	print "\n\nWITHDRAWN HARD CODED   Done\n ";
	return;
}

=head2 no_genotype_results_for_run_id

  Arg [1]   : DB adaptor
  Arg [2]   : verbose
 
  Function  : Should have genotype results for most things


=cut

sub no_genotype_results_for_run_id {
	my ( $db, $verbose ) = @_;

	my $sql =
	    'select run_id, sample_name, center_name  from run_meta_info ' .
            'where status = "public" and run_id not in (select name from genotype_results)';
#	print $sql;
	my $sth = $db->dbc->prepare($sql);
	$sth->execute;
	my $ctr = 0;
        my $no_snps;

	while ( my @row = $sth->fetchrow_array ) {
		
	 
	  print "@row" if( $verbose);
	  if (! defined $snps{$row[1]} ){
	    print " No snps available\n" if( $verbose);
	    $no_snps++;
	  }
	  else{
	    print " $row[1] ", $snps{$row[1]} ," \n" if( $verbose);
	  }
	 

		$ctr++;
	}
	$sth->finish;

	print "\nHave $ctr run_meta_info entries with no genotype results\n";
	print "Have $no_snps runs have no snps\n\n";	
	return;
}

=head2 more_read_than_bases

  Arg [1]   : DB adaptor
  Arg [2]   : verbose
 
  Function  : Occasionally read and base counts get flipped at ERA.


=cut

sub more_read_than_bases {

	my ( $db, $verbose ) = @_;

	my $sql =
	    "select "
	  . "run_id, sample_name, center_name, "
	  . "archive_read_count, archive_base_count "
	  . "from run_meta_info "
	  . "where archive_read_count > archive_base_count";

	my $sth = $db->dbc->prepare($sql);
	$sth->execute;
	my $ctr = 0;

	while ( my @row = $sth->fetchrow_array ) {
		print "@row\n" if $verbose;
		$ctr++;
	}
	$sth->finish;

	print "Have $ctr run_meta_info entries with read counts > base count\n";
	return;
}

=head2 fastq_attrib_stat_obj_check

  Arg [1]   : DB adaptor
  Arg [2]   : type
  Arg [3]   : verbose
 
  Function  : Make sure fastq files have base/read couhnts in Statistics table.


=cut

sub fastq_attrib_stat_obj_check {
	my ( $db, $type, $verbose ) = @_;

	my $sql = "select name, count(statistics.attribute_value) as cnt from file,
	statistics where file.file_id = statistics.other_id and
	statistics.table_name = 'file' and file.type = ? group
	by name having cnt = 1";

	my $sth = $db->dbc->prepare($sql);
	$sth->bind_param( 1, $type );
	$sth->execute;

	my $ctr = 0;

	while ( my @row = $sth->fetchrow_array ) {
		print "@row\n" if $verbose;
		$ctr++;
	}
	$sth->finish;

	print "Have $ctr fastq files of type $type";
	print " with only 1 attribute entry in Statistics table\n";
	return;
}

sub get_db_connection {
	my $o = shift;

	my $db = ReseqTrack::DBSQL::DBAdaptor->new(
		-host   => $$o{dbhost},
		-user   => $$o{dbuser},
		-port   => $$o{dbport},
		-dbname => $$o{dbname},
		-pass   => $$o{dbpass},
	);

	throw "No db connection" if ( !$db );

	#print "Got db connection\n";
	return $db;

}
