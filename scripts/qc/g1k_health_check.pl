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

GetOptions(
	\%input,    'dbhost=s', 'dbname=s',      'dbuser=s',
	'dbpass=s', 'dbport=s', 'working_dir=s', 'verbose!'
);


my $db = &get_db_connection( \%input );



print "\n\nStarting Jr Health check of: " . $input{dbname} . "\n";
print '-' x 45;
print "\n";


&history_objects_with_no_assoc_obj( $db, 'file');
&history_objects_with_no_assoc_obj( $db, 'collection');



&objects_with_no_history ( $db, 'collection');

&fastq_attrib_stat_obj_check( $db, 'FILTERED_FASTQ', $input{verbose} );

&fastq_attrib_stat_obj_check( $db,'WITHDRAWN_FILTERED_FASTQ',$input{verbose});


&sample_swap_results ( $db, $input{verbose});

&more_read_than_bases( $db, $input{verbose} );

&no_genotype_results_for_run_id( $db, $input{verbose} );

&check_fastq_collection_types_consistent( $db, 'FILTERED_FASTQ',
	$input{verbose} );

##&check_fastq_collection_types_consistent( $db, 'WITHDRAWN_FILTERED_FASTQ',
##	$input{verbose} );

print "\n\nDone \n";


=head2 history_objects_with_no_coll_obj 

  Arg [1]   : DB adaptor
  Arg [2]   : verbose
 
  Function  : collections get deleted, but not there history obj

=cut

sub history_objects_with_no_assoc_obj {
    my ( $db, $table,$verbose ) = @_;
    
    my $get = "select history_id,other_id, table_name from history where table_name = \"$table\" ";
    my $clause ="  other_id not in (select ${table}_id from $table)"; 

    my $sql = "$get AND  $clause" ;
    
    print $sql,"\n\n" if $verbose;
    
    my $sth = $db->dbc->prepare($sql);
    $sth->execute;
    my $ctr = 0;

    while ( my @row = $sth->fetchrow_array ) {
            $ctr++;
      }

 print "Got $ctr  history (table = $table) objects with no existing $table objects\n";
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
	    'select run_id, sample_name, center_name  from run_meta_info where status = "public" and run_id not in (select name from genotype_results)';
#	print $sql;
	my $sth = $db->dbc->prepare($sql);
	$sth->execute;
	my $ctr = 0;

	while ( my @row = $sth->fetchrow_array ) {
		print "@row\n";# if $verbose;
		$ctr++;
	}
	$sth->finish;

	print "Have $ctr run_meta_info entries with no genotype results\n";
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
