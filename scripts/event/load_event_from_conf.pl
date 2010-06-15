#!/sw/arch/bin/perl
use strict;
use ReseqTrack::Event;
use ReseqTrack::DBSQL::DBAdaptor;
use ReseqTrack::Tools::Exception;
use Getopt::Long;
use Data::Dumper;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;

my $file;
my $run;

my ($read, $write, $update, $load, $help, $verbose,);

&GetOptions(
    'dbhost=s'  => \$dbhost,
    'dbname=s'  => \$dbname,
    'dbuser=s'  => \$dbuser,
    'dbpass=s'  => \$dbpass,
    'dbport=s'  => \$dbport,
    'file=s'    => \$file,
    'read!'     => \$read,
    'write!'    => \$write,
    'help!'     => \$help,
    'run!'      => \$run,
    'verbose!'  => \$verbose,
);

if ($help) {
    perldocs();
}

#Could allow this. Add to db then dump contents. Always dump contents
#after modifying db.
if ($read && $write) {
    throw "Must specify -read or -write. Not both\n";
}

if (!$read && !$write) {
    print "Must specify -read or -write.\n";
    perldocs();
}

if (!$file) {
    throw "Must specify file (-file) to read from or write to";
}

my $db = ReseqTrack::DBSQL::DBAdaptor->new(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $dbname,
    -pass   => $dbpass,
);

if ($read) {
    my $events_list = parse_files([$file]);
    process_event_from_file($events_list, $verbose);
}

if ($write) {
    dump_event_table_contents($file, $verbose);
}

#######################################################
sub dump_event_table_contents {

    #output contents of 'event' table to fomatted text file

    my $outfile = shift;
    my $verbose = shift;

    my $ea     = $db->get_EventAdaptor;
    my $events = $ea->fetch_all;

    open(OUT, '>', $outfile) or die "File open error $outfile: $!";

    printf "Event data entries:%d \n", scalar @$events;

    foreach my $ievent (@$events) {
        printf OUT "[$ievent->{name}]\n"          if $ievent->{name};
        printf OUT "program=$ievent->{program}\n" if $ievent->{program};
        printf OUT "program_version=$ievent->{program_version}\n"
          if $ievent->{program_version};
        printf OUT "options=$ievent->{options}\n" if $ievent->{options};
        printf OUT "input_flag=$ievent->{input_flag}\n"
          if $ievent->{input_flag};
        printf OUT "farm_options=$ievent->{farm_options}\n"
          if $ievent->{farm_options};
        printf OUT "batch_size=$ievent->{batch_size}\n"
          if $ievent->{batch_size};
        printf OUT "output_path=$ievent->{output_path}\n"
          if $ievent->{output_path};
        printf OUT "type=$ievent->{type}\n" if $ievent->{type};
        printf OUT "table_name=$ievent->{table_name}\n"
          if $ievent->{table_name};
        printf OUT "\n";
    }

    close(OUT);
    print "\nFinished. See $outfile file\n\n";
}

sub process_event_from_file {

    #loop through array of hashes and load/update event info
    # into database

    my $events  = shift;     #array of hashes with event info
    my $verbose = shift;
    my $class   = "Event";
    my $updated = 0;
    my $stored  = 0;

    printf "\nFiles to process:%d \n", scalar @$events if $verbose;

    my $ea = $db->get_EventAdaptor;

    foreach my $jevent (@$events) {
      
      my $stored_event = $ea->fetch_by_name($jevent->name);

        if ($stored_event) {
            $jevent->dbID($stored_event->dbID);
            $jevent->adaptor($ea);
            $ea->update($jevent);
            $updated++;
        } else {
            $ea->store($jevent);
            $stored++;
        }
    }

    print "Events stored : $stored\n";
    print "Events updated: $updated\n";
}



sub parse_files {

  my $files = shift;

  my %headers;     # will store names of headers and number of keys for each
  
  my $hcounter = 0;
  my %horder; # stores order in which entries were read


  my $config = {};
  # read each file

  foreach my $file (@$files) {

    if (! -e $file) {
      throw("analysis file $file not found\n");
    }
    my $header = "";

    open (FILE, $file) or throw "Couldn't open file $file";
    while (<FILE>) {
      chomp();

      # Comment or blank line
      next if (/^\s$/ || /^\#/);

      # [HEADER]
      if (/^\[(.*)\]\s*$/) {         # $1 will be the header name, without the [] 
	$header = $1;
	$headers{$header} = 0;
        $horder{$header} = $hcounter++;
	#print "Reading stanza $header\n";
      } 

      # key=value
      if (/^([^=\s]+)\s*=\s*(.+?)\s*$/) {   # $1 = key, $2 = value

	my $key = lc($1);           # keys stored as all lowercase, values have case preserved
	my $value = $2;
	if (length($header) == 0) {
	  throw("Found key/value pair $key/$value outside stanza");
	}
	#print "Key: $key Value: $value\n"; 
      	
	# Check if this header/key is already defined
	if (exists($config->{$header}->{$key})) {
	  throw("$key is already defined for [$header]; cannot be redefined");
	} else {
	  # store them in the config hash
	  $config->{$header}->{$key} = $value;
	  #print "$header:$key=$value\n";
	  $headers{$header}++;  # will be used to check for headers with no keys
	}

      }

    } # while <FILE>

    close FILE;
  }

  my @analyses;
  # add a blank key/value for any headers that have no keys
  foreach my $h (sort { $horder{$a} <=> $horder{$b} }  keys (%headers)) {
    if (!$config->{$h}->{'type'}) {
      throw("you seem to have no type for $h can't ".
            "create a an event object without an type");
    }
    if (!$config->{$h}->{'output_path'}) {
      throw("you seem to have no output_path for $h can't ".
            "create a an event object without an output_dir");
    }
    if (!$config->{$h}->{'table_name'}) {
      throw("you seem to have no table name for $h can't ".
            "create a an event object without an table name");
    }
    
    my $analysis = ReseqTrack::Event->new
      (
       -program         => $config->{$h}->{program},
       -program_version => $config->{$h}->{program_version},
       -options      => $config->{$h}->{options},
       -name      => $h,
       -input_flag => $config->{$h}->{input_flag},
       -farm_options => $config->{$h}->{farm_options},
       -batch_size => $config->{$h}->{batch_size},
       -output_path => $config->{$h}->{output_path},
       -type => $config->{$h}->{type},
       -table_name => $config->{$h}->{table_name},
      );
    push(@analyses, $analysis);
  }

  return \@analyses;

}
sub perldocs {
    exec('perldoc', $0);
    exit(0);
}

=pod

=head1 NAME

ReseqTrack/scripts/event/load_event_from_conf.pl

=head1 SYNOPSIS

    This script takes a config file in the following format and parse it to create
    and load event objects into the event table

        [event_name]
        program=/path/to/program
        program_version=program version
        options=command line options for version
        input_flag=command line flag to associate with input
        farm_options=options for LSF submission of event
        batch_size=number of inputs to associate in a single LSF submission
        output_path=where to direct the stdout and stderr of the file
        type=the type of the input data required by the event
        table_name=the table name needed to retrieve the input data

    It can also dump existing events in this format


=head2 OPTIONS

        -dbhost, the name of the mysql-host
        -dbname, the name of the mysql database
        -dbuser, the name of the mysql user
        -dbpass, the database password if appropriate
        -dbport, the port the mysql instance is running on, this defaults to 4197 the 
         standard port for mysql-g1kdcc.ebi.ac.uk
        -file, name of file to read in or write out
        -read, binary flag to indicate the file specified by -file should be parsed and 
         loaded into the database
        -write, binary flag to indicate the file specified should be written in the given
          format based on the workflow objects already in the database
        -help, binary flag to indicate the help should be printed


=head1 Example:


$DB_OPTS="-dbhost mysql-host -dbuser rw_user -dbpass **** -dbport 4197 -dbname my_database"


    This dumps the existing event objects into a file of the given name:
    perl ReseqTrack/scripts/event/load_event_from_conf.pl  $DB_OPTS -file /path/to/event.conf -write
   
    
    This will load the event objects specifed in event.conf:
    perl ReseqTrack/scripts/event/load_event_from_conf.pl  $DB_OPTS -file /path/to/event.conf -read
    The write process will overwrite any existing file with the given name

=head2 Other useful scripts

    ReseqTrack/workflow/load_workflow_from_conf.pl is another script which works in a
    very similar manner but for workflow objects

=cut

