#!/sw/arch/bin/perl -w

use strict;
use ReseqTrack::Tools::Exception;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::GeneralUtils;
use Getopt::Long;

$| = 1;

my $new_index_file;
my $old_index_file;
my $dbhost;
my $dbuser;
my $dbpass;
my $dbport = 4197;
my $dbname;
my $module;
my $sequence;
my $alignment;
my $sequence_module = 'ReseqTrack::Tools::Statistics::SequenceIndexStatistics';
my $alignment_module = 'ReseqTrack::Tools::Statistics::AlignmentIndexStatistics';
my $help = 0;
&GetOptions( 
	    'dbhost=s'       => \$dbhost,
	    'dbname=s'       => \$dbname,
	    'dbuser=s'       => \$dbuser,
	    'dbpass=s'       => \$dbpass,
	    'dbport=s'       => \$dbport,
            'new_summary_file=s'=> \$new_index_file,
	    'old_summary_file=s' => \$old_index_file,
	    'module=s' => \$module,
	    'sequence!' => \$sequence,
	    'alignment!' => \$alignment,
	    'help!' => \$help,
	   );

if($help){
  useage();
}

unless($module){
  if($sequence){
    $module = $sequence_module;
  }elsif($alignment){
    $module = $alignment_module;
  }else{
    throw("Seem to have no module specified, not sure how to generate stats using ".
	  $new_index_file." and ".$old_index_file);
  }
}

my $file = $module;
$file =~ s/::/\//g;
eval{
  require "$file.pm";
};
if($@){
  throw("Couldn't require $file $@");
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
  -host   => $dbhost,
  -user   => $dbuser,
  -port   => $dbport,
  -dbname => $dbname,
  -pass   => $dbpass,
    );


my $stats = $module->new(
			 -db => $db,
			 -new_index => $new_index_file,
			 -old_index => $old_index_file,
			);

$stats->make_stats;
$stats->print_stats;


=pod

=head1 NAME

ReseqTrack/scripts/process/generate_stats.pl

=head1 SYNOPSIS

This script can use the ReseqTrack::Tools::Statistics modules to generate summary 
statisitcs on the basis of sequence indexes or bas files. Other stats module could be
written. They must implement two methods, make_stats to generate the statistics and
print_stats to print out a csv file.

=head1 OPTIONS

Database options

This is the database the index files will be fetched from if none is specified on the
command line

-dbhost, the name of the mysql-host
-dbname, the name of the mysql database
-dbuser, the name of the mysql user
-dbpass, the database password if appropriate
-dbport, the port the mysql instance is running on, this defaults to 4197 the standard
port for mysql-g1kdcc.ebi.ac.uk

Other options

-old_summary_file, path to the older index or bas file.
-new_summary_file, path to newer index or bas file. Index files must be named in the
form YYYYMMDD.sequence.index and if no files are passed in then index files are 
fetched from the given database and the 2 most recent files are worked out on the basis of the specified date.
-module, this is the perl path to the module which should be used e.g ReseqTrack::Tools::Statistics::SequenceIndexStatistics
-sequence, this means you default to using the SequenceIndexStatistics module
-alignment, this means you default to using the AlignmentIndexStatistics module

=head1 Examples

perl generate_stats.pl -dbhost myhost -dbuser ro -dbname mydatabase -sequence

this will automatically find the two most recent sequence index files and 
produce stats for them

=cut
