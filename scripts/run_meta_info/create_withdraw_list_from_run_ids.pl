# Use run_ids (combined with "type") to create a withdraw list.
# File produced is tab separated list that can be used to withdraw files
# from database. Pulls file names via Collection table.

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
use File::Path;
use Data::Dumper;

my $host_name;
my $file;    # unpaired fastq file
my $run_id;
my $type;
my $dbhost = "mysql-g1kdcc.ebi.ac.uk";
my $dbname = "g1k_archive_staging_track";
my $dbuser = "g1krw";
my $dbport = 4197;
my $dbpass = "thousandgenomes";
my $verbose;
my $run;
my $total_to_pull = 0;
my $outfile       = "withdraw.lst";
my $help;
my @check_ftp_exists;
my @check_withdrawn_exists;
my $check_exists = 0;

&GetOptions(

	'file=s'    => \$file,
	'run_id=s'  => \$run_id,
	'type=s'    => \$type,
	'dbhost=s'  => \$dbhost,
	'dbname=s'  => \$dbname,
	'dbuser=s'  => \$dbuser,
	'dbpass=s'  => \$dbpass,
	'dbport=s'  => \$dbport,
	'outfile=s' => \$outfile,
	'verbose'   => \$verbose,
	'check'     => \$check_exists,
	'help'      => \$help,
);

if ($help) {
	useage();
}

warning "No -type specified. Pulling all matching run_id entries" if ( !$type );
sleep(3) if ( !$type );

my @run_ids;

push( @run_ids, $run_id ) if ($run_id);

if ($file) {
	open( FH, '<', $file ) || die "Could not open $file";
	while (<FH>) {
		chomp;
		push( @run_ids, $_ );
	}
	close(FH);
}
print "\nNumber of run_ids to pull : ", scalar(@run_ids), "\n";

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
	-host   => $dbhost,
	-user   => $dbuser,
	-port   => $dbport,
	-dbname => $dbname,
	-pass   => $dbpass,
);

open( FO, '>', $outfile ) || die "Could not open $outfile";

# query information from Collection table
my $ca = $db->get_CollectionAdaptor;
my $collection;

foreach my $i (@run_ids) {

	if ($type) {
		$collection = $ca->fetch_by_name_and_type( $i, $type );
	}
	else {
		$collection = $ca->fetch_by_name($i);    #collection array
	}

	if ( !defined($collection) ) {
		throw "Could not pull collection info for run_id = $i type= $type";
	}

	if ($type) {
		my $others = $collection->others;

		foreach my $q (@$others) {
			print $q->{name}, "\t $q->{type}\n" if $verbose;
			my $name = $q->{name};
			push( @check_ftp_exists, $name );
			print FO $name, "\t";
			$name =~ s/ftp/withdrawn/;
			push( @check_withdrawn_exists, $name );
			print FO $name, "\n";
			$total_to_pull++;
		}

	}
	else {
		foreach my $q (@$collection) {
			my $others = $q->others;
			foreach my $f (@$others) {
				print "$f->{name} \t $f->{type}\n" if $verbose;
				my $name = $f->{name};
				push( @check_ftp_exists, $name );
				print FO $name, "\t";
				$name =~ s/ftp/withdrawn/;
				push( @check_withdrawn_exists, $name );
				print FO $name, "\n";
				$total_to_pull++;
			}
		}
	}
}
close(FO);

print "Total number of files to pull = $total_to_pull\n";
print "See $outfile for move list\n\n";


#Just checking to see if all files are there on ftp or
#might get clobbered in withdrawn.'

if ($check_exists) {
	my $ctr = 0;
	print scalar(@check_ftp_exists), "\n" if $verbose;

	foreach my $i (@check_ftp_exists) {
		print "Does not exist: $i\n" if ( (!-e $i)  && $verbose );
		$ctr++;
	}
	if ( $ctr){
		print "Some file do not exist: Total $ctr\n\n";
		
	}
	
 	$ctr = 0;
	print scalar(@check_withdrawn_exists), "\n" if $verbose;
	foreach my $j (@check_withdrawn_exists) {
		if ( !-e $j ) {
			print "Does not exist: $j\n" if $verbose;
			$ctr++;
		}
	}
	
	if ( $ctr != scalar(@check_withdrawn_exists) ) {
		print "It appears that some files may get over written\n";
	}

}





sub useage {
	exec( 'perldoc', $0 );
	exit(0);
}

=pod

=head1 NAME

		ReseqTrack/scripts/run_meta_info/create_withdraw_list_from_run_ids.pl

=head1 SYNOPSIS
		Create a move list to withdraw files based on run_id (and type). You cannot
		just give type. It will not work.


=head1 OPTIONS
		
		-dbhost, the name of the mysql-host
		-dbname, the name of the mysql database
		-dbuser, the name of the mysql user
		-dbpass, the database password if appropriate
		-dbport, the port the mysql instance is running on, this defaults to 4197 the
 		standard port for mysql-g1kdcc.ebi.ac.uk
		-file, input name of run_ids you wish to withdraw
		-run_id,  sinlge run_id to withdraw
		-type type of file associated with run_id
		-help, Binary flag to indicate if the help should be printed out
		-verbose.


=head1 Examples
		
perl $RESEQTRACK/scripts/ create_withdraw_list_from_run_ids.pl 
-dbhost mysql-g1kdcc.ebi.ac.uk -dbport 4197 -dbuser g1krw -dbpass xxxxxxxxxx
-dbname g1k_archive_staging_track -file bad_run_ids.lst -type FILTERED_FASTQ
    
perl $RESEQTRACK/scripts/ create_withdraw_list_from_run_ids.pl 
-dbhost mysql-g1kdcc.ebi.ac.uk -dbport 4197 -dbuser g1krw -dbpass xxxxxxxxxx
-dbname g1k_archive_staging_track -file bad_run_ids.lst 
   
perl $RESEQTRACK/scripts/ create_withdraw_list_from_run_ids.pl 
-dbhost mysql-g1kdcc.ebi.ac.uk -dbport 4197 -dbuser g1krw -dbpass xxxxxxxxxx
-dbname g1k_archive_staging_track -run_id SRR013914 -type FILTERED_FASTQ
    
perl $RESEQTRACK/scripts/ create_withdraw_list_from_run_ids.pl 
-dbhost mysql-g1kdcc.ebi.ac.uk -dbport 4197 -dbuser g1krw -dbpass xxxxxxxxxx
-dbname g1k_archive_staging_track -run_id SRR013914 
    
=cut

