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
use ReseqTrack::DBSQL::HistoryAdaptor;
use ReseqTrack::History;
use File::Basename;
use File::Path;
use Data::Dumper;

my $host_name;
my $file;    # unpaired fastq file
my $list;
my $type;
my $dbhost;
my $dbname;
my $dbuser;
my $dbport;
my $dbpass;
my $help;
my $verbose;
my $new_type;;
my $comment;
my $run;

&GetOptions(
	'file=s'      => \$file,
	'file_list=s' => \$list,
	'type=s'      => \$type,
	'dbhost=s'    => \$dbhost,
	'dbname=s'    => \$dbname,
	'dbuser=s'    => \$dbuser,
	'dbpass=s'    => \$dbpass,
	'dbport=s'    => \$dbport,
	'comment=s'   => \$comment,
	'new_type=s'  => \$new_type,
	'help'        => \$help,
	'run' => \$run,
	'verbose'        => \$verbose,
);

throw "new type not set ( -new_type)" if (! $new_type);
throw "history comment  not set ( -comment)" if (! $comment);



my @files;

push( @files, $file ) if ($file);

if ($list) {
	open( FH, '<', $list ) || die "Could not open $list";
	while (<FH>) {
		chomp;
		push( @files, $_ );
	}
	close(FH);
}
print "\nNumber of file objects to pull : ", scalar(@files), "\n";
print "\nNOT IN RUN MODE. NOTHING WILL HAPPEN\n" if (! $run);

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
	-host   => $dbhost,
	-user   => $dbuser,
	-port   => $dbport,
	-dbname => $dbname,
	-pass   => $dbpass,
);

my $fa = $db->get_FileAdaptor;

foreach my $f (@files) {	
	print $f,"\n" if $verbose;
	my $i = $fa->fetch_by_name($f);
  	my $j = copy_file_object ($i, 1);
	$j->{type} = $new_type;
	my $new_history = create_history( undef, $i, $comment );
	$j->history($new_history);
	$fa->store( $j, 1 ) if $run;
}



